function qc = collect_isolation_qc(recording,signals,origins,pool,noisePool,ap,noise,result,d)
% Compact diagnostics from data already loaded/filtered. No RNG or rescoring.
% The displayed waveforms are a deterministic subset (at most 120 each);
% Waveform/PCA panels use bounded display subsets; affinities use every AP.
qc=d.QC; fs=recording.SampleRateHz; maxShown=120;
qc.TargetWaveforms=ap(display_rows(size(ap,1),maxShown),:);
qc.NoiseWaveforms=noise(display_rows(size(noise,1),maxShown),:);
qc.TargetMean=d.AllAPMean; qc.NoiseMean=mean(noise,1);
qc.TargetTimeSeconds=(pool(:,2)-1)/fs; qc.TargetAmplitude=zeros(size(pool,1),1);
qc.TargetSegment=pool(:,1); intervals=cell(size(d.Segments,1),1);
% Amplitude uses aligned troughs; ISIs use the original sorted timestamps.
% Only timestamps within the same stable/file intersection form intervals.
for s=1:numel(signals)
    rows=find(pool(:,1)==s);
    qc.TargetAmplitude(rows)=-reshape(signals{s}(pool(rows,2)-origins(s)),[],1);
    times=recording.SpikeSamples(recording.SpikeSamples>=d.Segments(s,1) & recording.SpikeSamples<=d.Segments(s,2));
    intervals{s}=diff(times(:))/fs*1000; % Never span a stable-block or file gap.
end
qc.ISIMilliseconds=vertcat(intervals{:});
qc.ShortISIReferenceMs=2;
qc.ShortISIFraction=mean(qc.ISIMilliseconds<qc.ShortISIReferenceMs);
shownAP=ap(display_rows(size(ap,1),3000),:);
shownNoise=noise(display_rows(size(noise,1),3000),:);
x=[shownAP;shownNoise]; centered=x-mean(x,1);
% Joint PCA is a visualization only. The metric uses the full waveform grid.
[~,singular,coeff]=svd(centered,'econ');
variance=diag(singular).^2;
if numel(variance)>=2 && sum(variance)>0
    scores=centered*coeff(:,1:2);
    qc.PCATarget=scores(1:size(shownAP,1),:); qc.PCANoise=scores(size(shownAP,1)+1:end,:);
    qc.PCAExplained=100*variance(1:2)/sum(variance);
else
    qc.PCATarget=[]; qc.PCANoise=[]; qc.PCAExplained=[NaN;NaN];
end
qc.TargetAffinity=[];
if isfield(d.Metric,'TargetAffinity'), qc.TargetAffinity=d.Metric.TargetAffinity; end
qc.Counts=[result.InputSpikes,result.EligibleSpikes,result.NumAPs,result.EligibleNoise,result.NumNoise];
% These checks link the plotted populations back to the numerical result.
events=[d.SelectedAP;d.SelectedNoise]; bounds=d.Segments(events(:,1),1:2);
qc.CrossBoundaryWindows=sum(events(:,2)-d.Config.BinSamples<bounds(:,1) | events(:,2)+d.Config.BinSamples>bounds(:,2));
qc.DisplayCounts=[size(ap,1),size(noise,1)];
qc.CountsMatch=size(pool,1)==result.NumAPs && size(noisePool,1)==result.NumNoise;
qc.DisplayCountsMatch=size(ap,1)==size(d.SelectedAP,1) && size(noise,1)==size(d.SelectedNoise,1);
qc.RepresentativeIS=mean(qc.TargetAffinity); qc.AggregateIS=result.IS;
qc.AffinityMeanError=abs(qc.RepresentativeIS-result.IS);
% One illustrative trace around the middle eligible target, not an exhaustive audit.
chosen=pool(ceil(size(pool,1)/2),:); s=chosen(1); center=chosen(2);
halfWidth=round(.025*fs);
first=max(d.Segments(s,1),center-halfWidth); last=min(d.Segments(s,2),center+halfWidth);
range=first:last; sig=signals{s};
targets=pool(pool(:,1)==s & pool(:,2)>=first & pool(:,2)<=last,2);
others=noisePool(noisePool(:,1)==s & noisePool(:,2)>=first & noisePool(:,2)<=last,2);
input=recording.SpikeSamples(recording.SpikeSamples>=first & recording.SpikeSamples<=last);
qc.Trace=struct('TimeMs',(range-first)/fs*1000,'Signal',sig(range-origins(s)), ...
    'TargetMs',(targets-first)/fs*1000,'TargetAmplitude',reshape(sig(targets-origins(s)),[],1), ...
    'NoiseMs',(others-first)/fs*1000,'NoiseAmplitude',reshape(sig(others-origins(s)),[],1), ...
    'InputMs',(input-first)/fs*1000,'Segment',s,'FileIndex',d.Segments(s,3), ...
    'GlobalSamples',[first last],'Threshold',d.Threshold);
end

function rows=display_rows(n,limit)
rows=unique(round(linspace(1,n,min(n,limit))));
end
