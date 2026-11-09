function events = threshold_event_centers(signal,threshold,maxGap)
% Use the first minimum in each threshold-excursion group to resolve ties.
% Gap samples are above threshold, so cannot contain the group's minimum.
crossings=find(signal(:)<threshold);
if isempty(crossings), events=zeros(0,1); return; end
groups=cumsum([1;diff(crossings)>maxGap]);
values=signal(crossings); values=values(:);
minima=accumarray(groups,values,[],@min);
isMinimum=values==minima(groups);
events=accumarray(groups(isMinimum),crossings(isMinimum),[],@min);
end
