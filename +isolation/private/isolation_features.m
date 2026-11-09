function out = isolation_features(signals,origins,events,bin,split)
% Shared AP/noise grid: spline interpolation, local trough, common crop,
% then per-waveform mean subtraction. Call in bounded event batches.
offsets=-round(bin/3)*split:(bin-round(bin/3))*split-1;
grid=-bin:1/split:bin; center=bin*split+1;
out=zeros(size(events,1),numel(offsets));
for s=unique(events(:,1))'
    rows=find(events(:,1)==s); ix=events(rows,2)-origins(s);
    snippets=signals{s}(ix+(-bin:bin));
    interpolated=interp1(-bin:bin,snippets',grid,'spline');
    up=reshape(interpolated,numel(grid),numel(rows))';
    [~,k]=min(up(:,center-split:center+split),[],2);
    centers=center-split+k-1;
    cols=centers+offsets;
    linear=sub2ind(size(up),repmat((1:numel(rows))',1,numel(offsets)),cols);
    values=up(linear); out(rows,:)=values-mean(values,2);
end
end
