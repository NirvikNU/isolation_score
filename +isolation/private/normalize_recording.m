function [internal,mask,id] = normalize_recording(recording)
% Validate the public contract before calling the unchanged recording core.
% Separate stable intervals must NEVER become one consecutive stable trial
% run. Insert explicitly unstable gap rows for the internal trial interface.
validateattributes(recording,{'struct'},{'scalar'});
required={'Signal','SpikeSamples','SampleRateHz','StableIntervals'};
allowed=[required,{'FileBounds','ID'}];
assert(all(isfield(recording,required)),'isolation:RecordingFields', ...
    'Recording needs Signal, SpikeSamples, SampleRateHz and StableIntervals.');
unknown=setdiff(fieldnames(recording),allowed);
assert(isempty(unknown),'isolation:RecordingFields','Unknown recording fields: %s.',strjoin(unknown,', '));
validateattributes(recording.Signal,{'numeric'},{'real','vector','nonempty','nonsparse'});
n=numel(recording.Signal);
validateattributes(recording.SampleRateHz,{'numeric'},{'real','finite','scalar','>',600});
spikes=recording.SpikeSamples;
validateattributes(spikes,{'numeric'},{'real','finite','integer','positive','nonsparse'});
assert(isempty(spikes) || isvector(spikes),'isolation:SpikeIndices','SpikeSamples must be a vector.');
spikes=double(spikes(:));
assert(all(spikes<=n) && all(diff(spikes)>0),'isolation:SpikeIndices', ...
    'SpikeSamples must be sorted, unique, one-based sample indices within Signal.');
stable=intervals(recording.StableIntervals,n,'StableIntervals');
files=[1 n];
if isfield(recording,'FileBounds'), files=intervals(recording.FileBounds,n,'FileBounds'); end
assert(~isempty(files) && files(1,1)==1 && files(end,2)==n && ...
    all(files(2:end,1)==files(1:end-1,2)+1),'isolation:FileCoverage', ...
    'FileBounds must partition Signal without overlaps or uncovered samples.');
% Adjacent retained intervals are physically continuous, so merge them.
merged=zeros(size(stable)); count=0;
for k=1:size(stable,1)
    if count>0 && stable(k,1)==merged(count,2)+1
        merged(count,2)=stable(k,2);
    else
        count=count+1; merged(count,:)=stable(k,:);
    end
end
stable=merged(1:count,:);
if count==0
    trial=[1 n]; mask=false;
else
    trial=zeros(2*count-1,2); mask=false(2*count-1,1);
    trial(1:2:end,:)=stable; mask(1:2:end)=true;
    if count>1, trial(2:2:end,:)=[stable(1:end-1,2)+1,stable(2:end,1)-1]; end
end
internal=struct('Raw',recording.Signal,'SpikeSamples',spikes, ...
    'SampleRateHz',recording.SampleRateHz,'TrialSamples',trial,'FileBounds',files);
id="unit";
if isfield(recording,'ID')
    value=recording.ID;
    assert((ischar(value) && isrow(value)) || (isstring(value) && isscalar(value)) || ...
        (isnumeric(value) && isscalar(value) && isfinite(value) && isreal(value)), ...
        'isolation:ID','ID must be a text scalar or finite numeric scalar.');
    id=string(value);
    assert(~ismissing(id) && strlength(id)>0,'isolation:ID','ID must not be empty.');
end
end

function bounds=intervals(value,n,name)
validateattributes(value,{'numeric'},{'real','finite','integer','positive','2d','ncols',2,'nonsparse'},mfilename,name);
bounds=double(value);
assert(all(bounds(:,2)>=bounds(:,1)) && all(bounds(:,2)<=n) && ...
    all(bounds(2:end,1)>bounds(1:end-1,2)),'isolation:Intervals', ...
    '%s must contain sorted, disjoint inclusive intervals within Signal.',name);
end
