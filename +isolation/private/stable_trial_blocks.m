function [blocks,segments,trialRuns] = stable_trial_blocks(trialBounds,stableMask,fileBounds)
% Inclusive one-based sample bounds, including intertrial data within runs.
% segments splits stable blocks at acquisition-file boundaries (column 3=file).
% trialBounds is trials x 2 (events 1 and 11); stableMask follows row order.
% trialRuns gives inclusive first/last trial rows. blocks gives their outer
% sample bounds, retaining intertrial samples within each consecutive run.
% segments is [first sample, last sample, file row]; never filter or compute
% cross-event intervals across different segment rows.
validateattributes(trialBounds,{'numeric'},{'2d','ncols',2,'real','finite','integer','positive'});
validateattributes(fileBounds,{'numeric'},{'2d','ncols',2,'real','finite','integer','positive'});
assert(isvector(stableMask) && numel(stableMask)==size(trialBounds,1), ...
    'isolation:TrialCount','Stable-mask length must equal trial count.');
assert(all(ismember(stableMask(:),[0 1])),'isolation:StableMask','Use binary indicators.');
assert(all(trialBounds(:,2)>=trialBounds(:,1)) && all(fileBounds(:,2)>=fileBounds(:,1)), ...
    'isolation:IntervalBounds','Every trial/file interval must end at or after its start.');
assert(all(trialBounds(2:end,1)>trialBounds(1:end-1,2)), ...
    'isolation:TrialOverlap','Trial sample intervals overlap or are out of order.');
assert(all(fileBounds(2:end,1)==fileBounds(1:end-1,2)+1), ...
    'isolation:FileBounds','File coordinates must form the concatenated recording.');
assert(all(trialBounds(:,1)>=fileBounds(1,1)) && all(trialBounds(:,2)<=fileBounds(end,2)), ...
    'isolation:TrialRange','All trial bounds must fall inside the recording file bounds.');
changes=diff([false;logical(stableMask(:));false]);
trialRuns=[find(changes==1),find(changes==-1)-1];
blocks=[trialBounds(trialRuns(:,1),1),trialBounds(trialRuns(:,2),2)];
segments=zeros(0,3);
for b=1:size(blocks,1)
    for f=1:size(fileBounds,1)
        first=max(blocks(b,1),fileBounds(f,1)); last=min(blocks(b,2),fileBounds(f,2));
        if first<=last, segments(end+1,:)=[first,last,f]; end %#ok<AGROW>
    end
end
end
