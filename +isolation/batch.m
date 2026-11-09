function [results,diagnostics] = batch(jobs,options)
%BATCH Score ordered unit jobs, optionally reusing complete saved results.
%   Each job has ID and exactly one of Recording or Load (zero-argument
%   function handle returning the public recording struct). Optional Result
%   and Diagnostics reuse prior work without loading the signal.
%   Outputs are equally ordered cell arrays. Save them using MATLAB save.
%   options.Score = scoring overrides; options.Progress = true/false;
%   options.Recompute = false; options.ErrorMode = 'continue' or 'error'.
%   CaptureQC is per job (default false), independent of summary plotting.
if nargin<2, options=struct; end
validateattributes(jobs,{'struct'},{'vector'});
validateattributes(options,{'struct'},{'scalar'});
cfg=struct('Score',struct,'Progress',true,'Recompute',false,'ErrorMode','continue');
names=fieldnames(options);
for k=1:numel(names)
    assert(isfield(cfg,names{k}),'isolation:BatchOption','Unknown batch option: %s.',names{k});
    cfg.(names{k})=options.(names{k});
end
validateattributes(cfg.Progress,{'logical','numeric'},{'scalar','binary'});
validateattributes(cfg.Recompute,{'logical','numeric'},{'scalar','binary'});
cfg.ErrorMode=validatestring(cfg.ErrorMode,{'continue','error'});
scoreOptions=cfg.Score;
if ~isfield(scoreOptions,'Progress'), scoreOptions.Progress=cfg.Progress; end
isolation.defaults(scoreOptions); % Catch options errors once, before processing.
assert(~isfield(scoreOptions,'CaptureQC'),'isolation:BatchQC', ...
    'Set jobs(k).CaptureQC for selected units instead of Score.CaptureQC.');
allowed={'ID','Recording','Load','Result','Diagnostics','CaptureQC'};
assert(isempty(setdiff(fieldnames(jobs),allowed)),'isolation:JobFields','Unknown job field.');
n=numel(jobs); results=cell(size(jobs)); diagnostics=cell(size(jobs));
reuse=false(size(jobs)); ids=strings(size(jobs));
for k=1:n
    assert(isfield(jobs,'ID'),'isolation:JobID','Each job needs an ID.');
    ids(k)=string(jobs(k).ID);
    assert(~ismissing(ids(k)) && strlength(ids(k))>0,'isolation:JobID','Use nonempty scalar IDs.');
    reuse(k)=~cfg.Recompute && isfield(jobs,'Result') && complete(jobs(k).Result);
    if isfield(jobs,'CaptureQC') && ~isempty(jobs(k).CaptureQC)
        validateattributes(jobs(k).CaptureQC,{'logical','numeric'},{'scalar','binary'});
    end
end
assert(numel(unique(ids))==n,'isolation:JobID','Use distinct IDs, e.g. session/unit.');
total=sum(~reuse); finished=0; seconds=0;
if cfg.Progress, fprintf('Isolation: %d to calculate; %d complete results retained.\n',total,sum(reuse)); end
for k=1:n
    if reuse(k)
        results{k}=jobs(k).Result;
        if isfield(jobs,'Diagnostics') && isstruct(jobs(k).Diagnostics) && isscalar(jobs(k).Diagnostics)
            diagnostics{k}=jobs(k).Diagnostics;
        else
            diagnostics{k}=struct('ID',ids(k),'Result',results{k});
        end
        diagnostics{k}.ID=ids(k); diagnostics{k}.Result=results{k};
        diagnostics{k}.ReusedExisting=true;
        continue;
    end
    if cfg.Progress
        eta='estimating after first unit';
        if finished>0, eta=sprintf('about %.1f min remaining',(total-finished)*seconds/finished/60); end
        fprintf('Calculating %s (%d/%d); %s\n',ids(k),finished+1,total,eta);
    end
    timer=tic;
    try
        hasRecording=isfield(jobs,'Recording') && ~isempty(jobs(k).Recording);
        hasLoader=isfield(jobs,'Load') && ~isempty(jobs(k).Load);
        assert(xor(hasRecording,hasLoader),'isolation:JobInput','Provide exactly one of Recording or Load.');
        if hasLoader
            assert(isa(jobs(k).Load,'function_handle'),'isolation:JobInput','Load must be a function handle.');
            recording=jobs(k).Load();
        else
            recording=jobs(k).Recording;
        end
        recording.ID=ids(k); unitOptions=scoreOptions;
        unitOptions.CaptureQC=isfield(jobs,'CaptureQC') && ~isempty(jobs(k).CaptureQC) && jobs(k).CaptureQC;
        [results{k},diagnostics{k}]=isolation.score(recording,unitOptions);
        diagnostics{k}.ReusedExisting=false;
    catch error
        if strcmp(cfg.ErrorMode,'error'), rethrow(error); end
        results{k}=struct('IS',NaN,'snrAP',NaN,'snrBeforeAp',NaN,'fnScore',NaN,'fpScore',NaN, ...
            'NumAPs',0,'NumNoise',0,'Status','error','MetricVersion','Joshua2007-full-v1');
        diagnostics{k}=struct('ID',ids(k),'Result',results{k},'ReusedExisting',false, ...
            'ErrorIdentifier',error.identifier,'ErrorMessage',error.message);
        warning('isolation:JobFailed','%s: %s',ids(k),error.message);
    end
    elapsed=toc(timer); finished=finished+1; seconds=seconds+elapsed;
    diagnostics{k}.ElapsedSeconds=elapsed;
    if cfg.Progress, fprintf('  %s: IS %.4f; %s (%.1f s)\n',ids(k),results{k}.IS,results{k}.Status,elapsed); end
end
end

function yes=complete(result)
% A resume check is structural, not scientific validation of an older result.
fields={'IS','snrAP','snrBeforeAp','fnScore','fpScore'};
yes=isstruct(result) && isscalar(result) && all(isfield(result,fields));
if ~yes, return; end
for k=1:numel(fields)
    value=result.(fields{k});
    if ~(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value)), yes=false; return; end
end
end
