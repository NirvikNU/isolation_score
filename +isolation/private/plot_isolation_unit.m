function fig = plot_isolation_unit(detail,visible)
% Per-neuron checks from a saved diagnostics cell (no raw recording reload).
% CaptureQC selects which raw-data results retain these diagnostics.
% Plotting does not modify the scoring population or the RNG.
if nargin<2, visible='on'; end
name=char(detail.ID);
fig=figure('Color','w','Visible',visible,'Name',name,'Position',[25 25 1620 1110]);
if ~isfield(detail,'Scoring') || ~isfield(detail.Scoring,'QC') || ...
        ~isfield(detail.Scoring.QC,'TargetWaveforms')
    axis off;
    message=sprintf('%s\nStatus: %s\nNo waveform QC is available for this result.',name,detail.Result.Status);
    if isfield(detail,'PlotNote'), message=sprintf('%s\n%s',message,detail.PlotNote); end
    if isfield(detail,'ErrorMessage'), message=sprintf('%s\n%s',message,detail.ErrorMessage); end
    text(.05,.8,message,'Interpreter','none','FontSize',14);
    return;
end
d=detail.Scoring; q=d.QC; r=detail.Result; t=d.WaveformTimeMs;
% Display the sign convention: positive-going source waveforms are inverted
% together with noise for scoring. Older saved diagnostics remain supported.
polarityText=sprintf('Sign multiplier %+d for AP and noise',d.PolarityMultiplier);
if isfield(d,'PolaritySelection')
    p=d.PolaritySelection;
    polarityText=sprintf('Source phase: %s (%s); scoring sign %+d',p.Selected,p.Requested,d.PolarityMultiplier);
end
teal=[0 .4 .45]; orange=[.8 .32 .12];
layout=tiledlayout(fig,3,3,'Padding','compact','TileSpacing','loose');
title(layout,sprintf('%s | IS %.3f | %s',name,r.IS,strrep(r.Status,'_',' ')),'Interpreter','none','FontSize',17);
ax=nexttile(layout); hold(ax,'on'); targetWaveformAxis=ax;
% First row: waveforms actually scored, plus their shared display projection.
plot(ax,t,q.TargetWaveforms','Color',[.78 .86 .87],'LineWidth',.4);
plot(ax,t,q.TargetMean,'Color',teal,'LineWidth',2.3); xline(ax,0,':');
title(ax,sprintf('Display subset: %d shown / %d APs',size(q.TargetWaveforms,1),r.NumAPs));
xlabel(ax,'Time relative to normalized trough (ms)'); ylabel(ax,'Polarity-normalized amplitude'); grid(ax,'on');
ax=nexttile(layout); hold(ax,'on');
plot(ax,t,q.NoiseWaveforms','Color',[.91 .81 .74],'LineWidth',.4);
plot(ax,t,q.NoiseMean,'Color',orange,'LineWidth',2.3); xline(ax,0,':');
title(ax,sprintf('Display subset: %d shown / %d noise',size(q.NoiseWaveforms,1),r.NumNoise));
xlabel(ax,'Time relative to normalized trough (ms)'); ylabel(ax,'Polarity-normalized amplitude'); grid(ax,'on');
% A shared range makes AP/noise amplitudes directly comparable; include the
% full automatic range of both panels, then keep them linked during zooming.
sharedAmplitudeLimits=[min(targetWaveformAxis.YLim(1),ax.YLim(1)), ...
    max(targetWaveformAxis.YLim(2),ax.YLim(2))];
ylim(targetWaveformAxis,sharedAmplitudeLimits); ylim(ax,sharedAmplitudeLimits);
linkaxes([targetWaveformAxis ax],'y');
ax=nexttile(layout); hold(ax,'on');
if ~isempty(q.PCATarget)
    hN=scatter(ax,q.PCANoise(:,1),q.PCANoise(:,2),7,orange,'filled','MarkerFaceAlpha',.22);
    hT=scatter(ax,q.PCATarget(:,1),q.PCATarget(:,2),14,teal,'filled','MarkerFaceAlpha',.8);
    legend(ax,[hT hN],{'Target','Noise'},'Location','northwest','FontSize',9);
    % Keep the centre readable without hiding outliers: show full range inset.
    allPC=[q.PCATarget;q.PCANoise]; ordered=sort(allPC,1); n=size(ordered,1);
    lo=ordered(max(1,round(.005*n)),:); hi=ordered(min(n,round(.995*n)),:);
    padding=max((hi-lo)*.12,eps); lo=lo-padding; hi=hi+padding;
    outside=sum(any(allPC<lo | allPC>hi,2));
    xlim(ax,[lo(1) hi(1)]); ylim(ax,[lo(2) hi(2)]);
    pcaAxis=ax;
else
    outside=0; pcaAxis=[];
end
xlabel(ax,sprintf('PC1 (%.1f%%)',q.PCAExplained(1))); ylabel(ax,sprintf('PC2 (%.1f%%)',q.PCAExplained(2)));
title(ax,{'Shared PCA: display only',sprintf('Centre view: %d outside; full range inset',outside)}); grid(ax,'on');
ax=nexttile(layout); hold(ax,'on'); set(ax,'Color',[.94 .94 .94]);
% Second row: drift, gap-aware spike intervals, and contributions to mean IS.
top=max([1;q.TargetAmplitude])*1.05; bottom=min([0;q.TargetAmplitude])*1.05;
for b=1:size(d.StableBlocks,1)
    x=(d.StableBlocks(b,:)-1)/q.SampleRateHz/60;
    patch(ax,[x fliplr(x)],[bottom bottom top top],[.85 .94 .89],'EdgeColor','none','HandleVisibility','off');
end
for f=2:size(q.FileBounds,1)
    xline(ax,(q.FileBounds(f,1)-1)/q.SampleRateHz/60,':','Color',[.65 .65 .65],'HandleVisibility','off');
end
scatter(ax,q.TargetTimeSeconds/60,q.TargetAmplitude,3,teal,'filled','MarkerFaceAlpha',.28);
xlim(ax,[0 q.TotalSamples/q.SampleRateHz/60]); ylim(ax,[bottom top]);
title(ax,'Amplitude over retained blocks (green)');
xlabel(ax,'Concatenated file time (min; file edges dotted)'); ylabel(ax,'Aligned trough amplitude (sign flipped)'); grid(ax,'on');
ax=nexttile(layout);
histogram(ax,q.ISIMilliseconds,0:.2:10,'FaceColor',teal,'EdgeColor','none');
xline(ax,q.ShortISIReferenceMs,'--','Color',[.7 .2 .2]);
title(ax,sprintf('Within-segment ISIs: <2 ms = %.2f%%',100*q.ShortISIFraction));
xlabel(ax,'ISI (ms); intervals across gaps excluded'); ylabel(ax,sprintf('Count (all ISIs: %d)',numel(q.ISIMilliseconds)));
xlim(ax,[0 10]); grid(ax,'on');
ax=nexttile(layout);
if ~isempty(q.TargetAffinity)
    histogram(ax,q.TargetAffinity,0:.05:1,'FaceColor',teal,'EdgeColor','none'); hold(ax,'on');
    xline(ax,r.IS,'--','Color',[.7 .2 .2],'LineWidth',1.5);
end
xlim(ax,[0 1]); title(ax,sprintf('AP affinities; reported IS %.3f',r.IS));
xlabel(ax,'Target-neighbour affinity (not ground truth)'); ylabel(ax,'Target waveforms'); grid(ax,'on');
ax=nexttile(layout); hold(ax,'on'); tr=q.Trace;
% Third row: one trace, eligibility/scoring counts, and numerical checks.
hL=plot(ax,tr.TimeMs,tr.Signal,'Color',[.2 .25 .28],'LineWidth',.6);
hT=scatter(ax,tr.TargetMs,tr.TargetAmplitude,22,teal,'filled');
hN=scatter(ax,tr.NoiseMs,tr.NoiseAmplitude,14,orange);
hH=yline(ax,tr.Threshold,'--','Color',[.6 .2 .2],'LineWidth',1);
for x=tr.InputMs(:)'
    xline(ax,x,':','Color',[.7 .7 .7],'HandleVisibility','off');
end
legend(ax,[hL hT hN hH],{'Filtered','Aligned target','Noise event','Threshold'},'FontSize',8,'Location','best');
title(ax,sprintf('Example trace: file %d, segment %d',tr.FileIndex,tr.Segment));
xlabel(ax,'Snippet time (ms); dotted = input spike times'); ylabel(ax,'Polarity-normalized amplitude'); grid(ax,'on');
ax=nexttile(layout); bars=bar(ax,q.Counts,'FaceColor','flat');
bars.CData=[teal;teal;teal;orange;orange]; set(ax,'YScale','log');
xticks(ax,1:5); xticklabels(ax,{'In-block AP','Eligible AP','Scored AP','Eligible noise','Scored noise'}); xtickangle(ax,18);
for k=1:5
    if q.Counts(k)>0, text(ax,k,q.Counts(k)*1.18,sprintf('%d',q.Counts(k)),'HorizontalAlignment','center','FontSize',9); end
end
ylim(ax,[max(.5,min(q.Counts(q.Counts>0))*.5),max(q.Counts)*3]);
title(ax,'Events used for saved score'); ylabel(ax,'Event count (log scale)'); grid(ax,'on');
ax=nexttile(layout); axis(ax,'off');
lines={sprintf('SNR: waveform %.2f | pre-spike %.2f',r.snrAP,r.snrBeforeAp), ...
    sprintf('FN %.3f | FP %.3f',r.fnScore,r.fpScore), ...
    sprintf('Retained samples: %d',r.InputSamples), ...
    sprintf('%d biological blocks | %d file segments',size(d.StableBlocks,1),size(d.Segments,1)), ...
    sprintf('Alignment shift %d samples',d.GlobalShiftSamples),polarityText, ...
    sprintf('Waveform windows crossing boundaries: %d',q.CrossBoundaryWindows), ...
    sprintf('Counts match: %d | mean-affinity error %.2g',q.CountsMatch,q.AffinityMeanError)};
if ~isempty(d.Warnings), lines=[lines,d.Warnings]; end
if isfield(r,'MetricVersion')
    lines=[{r.MetricVersion},lines];
    lines=[{sprintf('Paper FN/FP; all-AP SNR uses %d APs',r.NumAPsSNR)},lines];
else
    lines=[{'Saved result has no metric version'},lines];
end
lines{end+1}='QC views support checks; they do not prove isolation.';
wrapped={};
for k=1:numel(lines)
    line=lines{k};
    while numel(line)>62
        space=find(line(1:62)==' ',1,'last'); if isempty(space), space=62; end
        wrapped{end+1}=line(1:space); line=strtrim(line(space+1:end)); %#ok<AGROW>
    end
    wrapped{end+1}=line; %#ok<AGROW>
end
text(ax,0,1,wrapped,'VerticalAlignment','top','Interpreter','none','FontSize',10);
title(ax,'Settings and numerical checks','HorizontalAlignment','left');
set(findall(fig,'Type','axes'),'FontSize',10,'Box','on');
if ~isempty(pcaAxis)
    drawnow;
    pos=pcaAxis.Position;
    inset=axes('Parent',fig,'Position',[pos(1)+.70*pos(3),pos(2)+.68*pos(4),.27*pos(3),.27*pos(4)]);
    scatter(inset,q.PCANoise(:,1),q.PCANoise(:,2),2,orange,'filled'); hold(inset,'on');
    scatter(inset,q.PCATarget(:,1),q.PCATarget(:,2),3,teal,'filled');
    set(inset,'FontSize',7,'Box','on'); title(inset,'Full range','FontSize',8);
end
end
