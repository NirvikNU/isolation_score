function fig = plotSummary(results,caption,visible)
%PLOTSUMMARY Summarize exactly the supplied results; never loads raw data.
%   Supply a cell array returned by isolation.batch, or {oneResult}.
%   Select cells to control sessions/groups. Unit plots are independent.
if nargin<2, caption='Isolation summary'; end
if nargin<3, visible='on'; end
if isstruct(results), results=num2cell(results); end
validateattributes(results,{'cell'},{'vector','nonempty'});
visible=validatestring(visible,{'on','off'});
fields={'IS','snrAP','snrBeforeAp','fnScore','fpScore'};
n=numel(results); values=nan(n,5); counts=nan(n,1); statuses=strings(n,1);
for k=1:n
    assert(isstruct(results{k}) && all(isfield(results{k},fields)), ...
        'isolation:Results','Every result must contain the five metric fields.');
    for j=1:5, values(k,j)=results{k}.(fields{j}); end
    if isfield(results{k},'NumAPs'), counts(k)=results{k}.NumAPs; end
    statuses(k)="saved_result";
    if isfield(results{k},'Status'), statuses(k)=string(results{k}.Status); end
end
fig=figure('Color','w','Visible',visible,'Name',char(caption),'Position',[50 50 1350 1000]);
layout=tiledlayout(fig,3,3,'Padding','compact','TileSpacing','loose');
title(layout,sprintf('%s | %d unit entries',caption,n),'Interpreter','none');
labels={'Isolation score','Waveform SNR','Pre-spike SNR','FN score','FP score'};
for j=1:5
    ax=nexttile(layout); finite=values(isfinite(values(:,j)),j);
    if ismember(j,[1 4 5]), histogram(ax,finite,0:.05:1);
    elseif ~isempty(finite), histogram(ax,finite,min(20,numel(finite))); end
    title(ax,sprintf('%s: finite N=%d',labels{j},numel(finite)));
    xlabel(ax,labels{j}); ylabel(ax,'Units'); grid(ax,'on');
end
ax=nexttile(layout); good=all(isfinite(values(:,1:2)),2);
scatter(ax,values(good,1),values(good,2),25,'filled');
xlabel(ax,'Isolation score'); ylabel(ax,'Waveform SNR'); grid(ax,'on');
ax=nexttile(layout); good=isfinite(values(:,1)) & isfinite(counts) & counts>0;
scatter(ax,counts(good),values(good,1),25,'filled'); set(ax,'XScale','log');
xlabel(ax,'AP count'); ylabel(ax,'Isolation score'); grid(ax,'on');
ax=nexttile(layout); [names,~,index]=unique(statuses,'stable');
barh(ax,accumarray(index,1)); yticks(ax,1:numel(names));
yticklabels(ax,strrep(names,'_',' ')); xlabel(ax,'Units'); title(ax,'Calculation status');
ax=nexttile(layout); axis(ax,'off');
text(ax,0,1,{'All supplied results are represented.', ...
    'NaN metrics are omitted from histograms.', ...
    'Inspect waveform and stability diagnostics.', ...
    'Low-population flags need inspection.', ...
    'FN/FP are classifier-based estimates.', ...
    'QC does not establish ground truth.'},'VerticalAlignment','top');
end
