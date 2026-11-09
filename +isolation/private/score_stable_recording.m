function [result,ap,noise,diagnostics] = score_stable_recording(recording,stableMask,options,replay)
% Score stable trial runs without filtering across trial-block or file joins.
% recording: Raw (vector), SpikeSamples (vector), TrialSamples (trials x 2),
% FileBounds (files x 2), SampleRateHz. Bounds are inclusive, one-based sample
% coordinates in Raw; stableMask follows TrialSamples row order.
% ap/noise have one waveform per row on the identical time grid. diagnostics
% retains original sample coordinates so a questionable event can be located.
if nargin<3, options=struct; end
if nargin<4, replay=[]; end
options.SampleRateHz=recording.SampleRateHz; cfg=isolation_config(options);
raw=recording.Raw(:)'; spikes=double(recording.SpikeSamples(:)');
validateattributes(spikes,{'numeric'},{'real','finite','integer','positive'});
assert(all(spikes<=numel(raw)) && all(diff(spikes)>0),'isolation:SpikeIndices','Spike indices must be unique, sorted, and in range.');
assert(recording.FileBounds(1,1)==1 && recording.FileBounds(end,2)==numel(raw), ...
    'isolation:FileCoverage','File bounds must cover Raw from sample 1 to its last sample.');
[blocks,segments,runs]=stable_trial_blocks(recording.TrialSamples,stableMask,recording.FileBounds);
result=struct('IS',NaN,'snrAP',NaN,'snrBeforeAp',NaN,'fnScore',NaN,'fpScore',NaN, ...
    'NumAPs',0,'NumNoise',0,'InputSamples',sum(blocks(:,2)-blocks(:,1)+1), ...
    'InputSpikes',0,'EligibleSpikes',0,'EligibleNoise',0,'Status','no_stable_trials', ...
    'NumAPsSNR',0,'NumBaselineSamples',0,'MetricVersion','Joshua2007-full-v1', ...
    'MetricStatistics',struct,'NoiseThreshold',NaN,'ThresholdFraction',cfg.ThresholdFraction);
ap=[]; noise=[];
diagnostics=struct('Config',cfg,'TrialRuns',runs,'StableBlocks',blocks,'Segments',segments, ...
    'Threshold',NaN,'GlobalShiftSamples',NaN,'PolarityMultiplier',NaN,'Warnings',{{}});
if isfield(recording,'TimestampQC') && recording.TimestampQC.ExactDuplicatesRemoved>0
    diagnostics.Warnings{end+1}=sprintf('Counted %d exactly repeated timestamps once across %d source files; see SourceFiles for times.', ...
        recording.TimestampQC.ExactDuplicatesRemoved,numel(recording.TimestampQC.AffectedExtensions));
end
if cfg.CaptureQC
    diagnostics.QC=struct('TrialSamples',recording.TrialSamples,'StableMask',logical(stableMask(:)), ...
        'FileBounds',recording.FileBounds,'TotalSamples',numel(raw),'SampleRateHz',recording.SampleRateHz);
end
if isempty(blocks), return; end
for b=1:size(blocks,1), result.InputSpikes=result.InputSpikes+sum(spikes>=blocks(b,1) & spikes<=blocks(b,2)); end
% Neither preprocessing nor full-population scoring changes the caller RNG.
bin=cfg.BinSamples; split=cfg.SplitRate; fs=cfg.SampleRateHz;
diagnostics.WaveformTimeMs=(-round(bin/3)*split:(bin-round(bin/3))*split-1)/split/fs*1000;
lagMax=round(cfg.ShiftSearchMs*fs/1000); alignMax=round(cfg.AlignmentSearchMs*fs/1000);
pad=ceil(cfg.FilterPaddingMs*fs/1000); guard=ceil(cfg.FileEdgeGuardMs*fs/1000);
margin=bin+lagMax+alignMax+1;
[order,cutoff]=buttord(300/(fs/2),100/(fs/2),3,10);
[bFilter,aFilter]=butter(order,cutoff,'high');
signals=cell(size(segments,1),1); origins=zeros(size(segments,1),1);
limits=zeros(size(segments,1),2); pool=zeros(0,2); known=cell(size(segments,1),1);
% Filter each stable/file intersection separately. Padding supplies real
% context from that file, but event eligibility stays inside the stable run.
% origins translates local filtered indices back to the original Raw vector.
% pool rows are [segment index, original sample]; limits include edge guards
% plus room for the full waveform and both alignment searches.
for s=1:size(segments,1)
    file=recording.FileBounds(segments(s,3),:);
    a=max(file(1),segments(s,1)-pad); z=min(file(2),segments(s,2)+pad);
    chunk=double(raw(a:z));
    assert(all(isfinite(chunk)),'isolation:NonfiniteSignal','Signal contains NaN or Inf in a selected interval.');
    if numel(chunk)<=3*max(numel(aFilter),numel(bFilter)), continue; end
    signals{s}=filtfilt(bFilter,aFilter,chunk); origins(s)=a-1;
    limits(s,:)=[max(segments(s,1),file(1)+guard)+margin,min(segments(s,2),file(2)-guard)-margin];
    selected=spikes(spikes>=limits(s,1) & spikes<=limits(s,2));
    pool=[pool;repmat(s,numel(selected),1),selected(:)]; %#ok<AGROW>
    known{s}=spikes(spikes>=segments(s,1) & spikes<=segments(s,2));
end
% Preserve the eligible count even when threshold estimation fails later.
result.EligibleSpikes=size(pool,1);
if size(pool,1)<2, result.Status='insufficient_eligible_spikes'; return; end
if cfg.Progress, fprintf('  Preparing all %d APs in %d stable/file segments...\n',size(pool,1),size(segments,1)); end
% All eligible APs determine preprocessing, independently of scoring budgets.
templateEvents=pool; ns=size(pool,1);
template=zeros(1,2*lagMax+1);
for s=1:numel(signals)
    ix=templateEvents(templateEvents(:,1)==s,2)-origins(s);
    if isempty(ix), continue; end
    for first=1:cfg.FeatureBatchSize:numel(ix)
        rows=first:min(first+cfg.FeatureBatchSize-1,numel(ix));
        template=template+sum(signals{s}(ix(rows)+(-lagMax:lagMax)),1);
    end
end
template=template/ns;
diagnostics.InitialTemplate=template;
positivePeak=max(template); negativeTrough=min(template);
if max(abs(template))==0, result.Status='flat_template'; return; end
% Choose one phase per unit/stable period, before thresholding or computing
% IS. Do not try both scores and keep the better one, or flip each waveform
% independently: either would change the populations being compared.
switch cfg.Polarity
    case 'auto', [~,peak]=max(abs(template)); polarity=-sign(template(peak));
    case 'negative', [~,peak]=min(template); polarity=1;
    case 'positive', [~,peak]=max(template); polarity=-1;
end
if polarity==0, result.Status='flat_template'; return; end
shift=peak-lagMax-1;
diagnostics.GlobalShiftSamples=shift; diagnostics.PolarityMultiplier=polarity;
selectedPolarity='negative'; if polarity==-1, selectedPolarity='positive'; end
peakMagnitudes=[max(0,positivePeak),max(0,-negativeTrough)];
balance=min(peakMagnitudes)/max(peakMagnitudes);
% This 90% flag is a descriptive QC heuristic, not a rejection criterion or
% a physiological waveform classification. Biphasic waveforms are allowed.
diagnostics.PolaritySelection=struct('Requested',cfg.Polarity,'Selected',selectedPolarity, ...
    'PositivePeak',positivePeak,'NegativeTrough',negativeTrough, ...
    'OppositeToDominantRatio',balance,'NearBalanced',balance>=0.9, ...
    'TemplateCount',ns,'PeakOffsetSamples',shift);
if strcmp(cfg.Polarity,'auto') && diagnostics.PolaritySelection.NearBalanced
    diagnostics.Warnings{end+1}=sprintf('Similar positive/negative template peaks; auto chose %s. Inspect biphasic waveform.',selectedPolarity);
end
% Normalize the ENTIRE filtered signal with this single sign. Target and
% noise extraction therefore share polarity, alignment and threshold sign.
for s=1:numel(signals), signals{s}=signals{s}*polarity; end
% Both populations use the same integer and sub-sample trough alignment.
for s=1:numel(signals)
    rows=find(pool(:,1)==s);
    if isempty(rows), continue; end
    pool(rows,2)=align_centers(signals{s},pool(rows,2)+shift-origins(s),alignMax)+origins(s);
end
% Stream all APs for SNR and threshold. Cache when within the memory budget.
% A memory-limited unit still receives its all-AP SNR and explicit IS status.
width=bin*split; cacheAP=8*ns*width<=cfg.MaxMemoryMB*1024^2/2;
apAll=[]; if cacheAP, apAll=zeros(ns,width); end
meanWave=zeros(1,width); m2=zeros(1,width); seen=0;
for first=1:cfg.FeatureBatchSize:ns
    rows=first:min(first+cfg.FeatureBatchSize-1,ns);
    batch=isolation_features(signals,origins,pool(rows,:),bin,split);
    if cacheAP, apAll(rows,:)=batch; end %#ok<AGROW> preallocated ns-by-width
    n=numel(rows); batchMean=mean(batch,1); difference=batchMean-meanWave;
    m2=m2+sum((batch-batchMean).^2,1)+difference.^2*(seen*n/(seen+n));
    meanWave=meanWave+difference*(n/(seen+n)); seen=seen+n;
end
result.NumAPsSNR=ns; p2p=max(meanWave)-min(meanWave);
noiseSD=sqrt(sum(m2)/(ns*width-1));
if noiseSD>0, result.snrAP=p2p/noiseSD; end
diagnostics.AllAPMean=meanWave;
[~,troughColumn]=min(meanWave);
if cacheAP
    values=apAll(:,troughColumn);
else
    values=zeros(ns,1);
    for first=1:cfg.FeatureBatchSize:ns
        rows=first:min(first+cfg.FeatureBatchSize-1,ns);
        batch=isolation_features(signals,origins,pool(rows,:),bin,split);
        values(rows)=batch(:,troughColumn);
    end
end
values=sort(values); weakCount=max(1,ceil(ns*cfg.ThresholdFraction));
threshold=mean(values(end-weakCount+1:end))/cfg.ThresholdRatio;
result.NoiseThreshold=threshold;
diagnostics.Threshold=threshold;
diagnostics.ThresholdSelection=struct('TemplateCount',ns,'WeakCount',weakCount, ...
    'TroughColumn',troughColumn,'NonnegativeCount',sum(values>=0), ...
    'WeakMean',threshold*cfg.ThresholdRatio,'UsesAllEligibleAPs',true);
if weakCount<20
    diagnostics.Warnings{end+1}=sprintf('Threshold weak tail has only %d APs; all eligible APs were used.',weakCount);
end
thresholdValid=isfinite(threshold) && threshold<0;
if ~thresholdValid
    result.Status='invalid_noise_threshold';
    diagnostics.Warnings{end+1}=sprintf(['Noise threshold %.6g is not negative: %d/%d eligible targets are nonnegative at the common trough. ' ...
        'Inspect waveform polarity/alignment before changing threshold settings.'],threshold,sum(values>=0),numel(values));
    if cfg.CaptureQC
        % Bounded failure evidence; no noise population or valid score exists.
        diagnostics.QC.ThresholdTargetMean=meanWave;
        if cacheAP
            [~,weakOrder]=sort(apAll(:,troughColumn),'descend');
            shown=weakOrder(1:min(120,weakCount));
            diagnostics.QC.ThresholdWeakWaveforms=apAll(shown,:);
            diagnostics.QC.ThresholdWeakEvents=pool(shown,:);
        end
    end
end
noisePool=zeros(0,2); baselineCount=0; baselineSum=0; baselineSumSq=0;
for s=1:numel(signals)
    if isempty(signals{s}) || limits(s,2)<limits(s,1), continue; end
    sig=signals{s};
    targetCenters=known{s}+shift-origins(s);
    targetCenters=targetCenters(targetCenters>alignMax & targetCenters<=numel(sig)-alignMax);
    targetCenters=align_centers(sig,targetCenters(:),alignMax);
    events=zeros(0,1);
    if thresholdValid, events=threshold_event_centers(sig,threshold,max(1,round(4*fs/30000))); end
    if ~isempty(events)
        globalEvents=events+origins(s);
        events=events(globalEvents>=limits(s,1) & globalEvents<=limits(s,2));
        tolerance=max(1,round(2*fs/30000));
        excluded=reshape(targetCenters+(-tolerance:tolerance),[],1);
        events=events(~ismember(events,excluded));
    end
    noisePool=[noisePool;repmat(s,numel(events),1),events+origins(s)]; %#ok<AGROW>
    % Exclude target spikes from the noise population with a small timing
    % tolerance. Other units' threshold events intentionally remain noise.
    % The pre-spike baseline is one BinSamples-wide window before each
    % isolated target; it cannot cross a stable/file boundary.
    targets=sort(targetCenters(:));
    isolated=targets([false;diff(targets)>3*bin]);
    g=isolated+origins(s);
    isolated=isolated(g-2*bin+1>=limits(s,1) & g-bin<=limits(s,2));
    if ~isempty(isolated)
        vals=sig(isolated+(-2*bin+1:-bin));
        baselineCount=baselineCount+numel(vals);
        baselineSum=baselineSum+sum(vals,'all'); baselineSumSq=baselineSumSq+sum(vals.^2,'all');
    end
end
result.EligibleNoise=size(noisePool,1);
if cfg.Progress, fprintf('  Populations ready: %d AP / %d noise\n',ns,result.EligibleNoise); end
result.NumBaselineSamples=baselineCount;
if baselineCount>1
    baselineSD=sqrt(max(0,(baselineSumSq-baselineSum^2/baselineCount)/(baselineCount-1)));
    if baselineSD>0, result.snrBeforeAp=p2p/baselineSD; end
end
if ~thresholdValid, return; end % SNR above does not depend on a noise threshold.
if ~cacheAP
    result.Status='waveform_memory_limit';
    diagnostics.Warnings{end+1}='All-AP SNR completed; AP waveform cache exceeds MaxMemoryMB/2.';
    return;
end
extractNoise=@(indices)isolation_features(signals,origins,noisePool(indices,:),bin,split);
if isempty(replay)
    [statistics,metric,ap,noise,apRows,noiseRows]=isolation_full_scores(apAll,size(noisePool,1),extractNoise,cfg);
else
    % Plot reconstruction repeats preprocessing only, never distance metrics.
    prior=replay.Scoring;
    assert(strcmp(replay.Result.MetricVersion,result.MetricVersion) && ...
        isequaln(prior.Threshold,threshold) && replay.Result.EligibleSpikes==ns && ...
        replay.Result.EligibleNoise==result.EligibleNoise,'isolation:PlotReplayMismatch', ...
        'Preprocessing/populations differ from the saved result.');
    statistics=replay.Result.MetricStatistics; metric=prior.Metric;
    apRows=prior.SelectedAPIndices; noiseRows=prior.SelectedNoiseIndices;
    ap=apAll(apRows,:); noise=extractNoise(noiseRows);
end
result.MetricStatistics=statistics;
names={'IS','fnScore','fpScore'};
for j=1:3, result.(names{j})=statistics.(names{j}).Mean; end
result.NumAPs=ns; result.NumNoise=size(noisePool,1);
diagnostics.Metric=metric;
diagnostics.SelectedAPIndices=apRows; diagnostics.SelectedNoiseIndices=noiseRows;
diagnostics.SelectedAP=pool(apRows,:); diagnostics.SelectedNoise=noisePool(noiseRows,:);
result.Status=metric.status;
if strcmp(result.Status,'ok') && min(ns,result.NumNoise)<cfg.MinimumReportedSpikes
    result.Status='ok_low_population';
end
for j=1:3
    st=statistics.(names{j});
    if ~isfinite(st.Mean) || st.LowPopulation
        diagnostics.Warnings{end+1}=sprintf('%s: %s (%d AP / %d noise).',names{j},st.Status,st.NumAPs,st.NumNoise);
    end
end
if cfg.CaptureQC && ~isempty(ap) && ~isempty(noise)
    diagnostics.QC=collect_isolation_qc(recording,signals,origins,pool,noisePool,ap,noise,result,diagnostics);
end
end

function centers=align_centers(signal,centers,radius)
% Integer trough alignment; the caller supplies centers with valid margins.
centers=centers(:);
if isempty(centers), return; end
windows=reshape(signal(centers+(-radius:radius)),numel(centers),2*radius+1);
[~,k]=min(windows,[],2);
centers=centers+k-radius-1;
end
