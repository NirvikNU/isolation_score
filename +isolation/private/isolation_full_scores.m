function [statistics,metric,ap,noise,apRows,noiseRows] = isolation_full_scores(apAll,noiseCount,extractNoise,cfg)
% Score every event once. Returned waveforms are bounded display subsets.
ns=size(apAll,1); width=size(apAll,2); nn=noiseCount;
apRows=display_rows(ns,cfg.DisplayEvents); noiseRows=display_rows(nn,cfg.DisplayEvents);
ap=apAll(apRows,:); noise=extractNoise(noiseRows);
pairs=double(ns)^2+double(ns)*nn+double(nn)^2;
native=strcmp(cfg.Backend,'native') || (strcmp(cfg.Backend,'auto') && pairs>cfg.MaxMatlabPairs);
if ns<2 || nn<1, native=false; end
if native
    assert(ispc && width==144 && cfg.K==31 && ns+nn>31, ...
        'isolation:NativeShape','Native requires Windows, 144 features, K=31 and >31 total events. Use MATLAB for small populations.');
    [fn,fp,is,metric]=native_scores(apAll,nn,extractNoise,cfg);
else
    assert(pairs<=cfg.MaxMatlabPairs || ns<2 || nn<1,'isolation:MatlabBudget', ...
        'Full MATLAB computation exceeds MaxMatlabPairs. Use native, or explicitly increase this guard.');
    assert(8*(numel(apAll)+nn*width)<cfg.MaxMemoryMB*1024^2/2, ...
        'isolation:MemoryBudget','MATLAB waveform arrays exceed the scoring memory budget.');
    allNoise=extractNoise(1:nn);
    if cfg.Progress, fprintf('  Full MATLAB metrics: %d AP / %d noise\n',ns,nn); end
    [fn,fp,is,metric]=knn_test(apAll,allNoise,cfg);
    metric.Backend='matlab'; metric.ISPrecision='FP64';
end
names={'IS','fnScore','fpScore'}; values=[is fn fp]; statistics=struct;
for j=1:3
    status='full_population'; if ~isfinite(values(j)), status=metric.status; end
    statistics.(names{j})=struct('Mean',values(j),'NumAPs',ns,'NumNoise',nn, ...
        'FullPopulation',true,'LowPopulation',min(ns,nn)<cfg.MinimumReportedSpikes, ...
        'Status',status,'Backend',metric.Backend,'Precision','FP64');
end
statistics.IS.Precision=metric.ISPrecision;
end

function [fn,fp,is,metric]=native_scores(a,nn,extractNoise,cfg)
ns=size(a,1); width=size(a,2); bytes=8*(numel(a)+nn*width);
zeroSpread=all(max(a,[],1)==min(a,[],1));
% Conservative scoring workspace estimate; raw channel memory is additional.
assert(3*bytes<cfg.MaxMemoryMB*1024^2,'isolation:MemoryBudget', ...
    'Estimated native workspace %.1f GB exceeds MaxMemoryMB.',3*bytes/1024^3);
nativeDir=fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'native');
assert(isfile(fullfile(nativeDir,'isolation_knn.exe')) && isfile(fullfile(nativeDir,'isolation_is.exe')), ...
    'isolation:NativeMissing','Native executables missing. See native/build.cmd.');
cuda=cfg.CUDARoot;
assert(zeroSpread || isfolder(fullfile(cuda,'bin')),'isolation:CUDA','Set options.CUDARoot or CUDA_PATH to the CUDA 11 toolkit folder.');
folder=tempname(cfg.TemporaryDirectory); mkdir(folder);
clean=onCleanup(@()rmdir(folder,'s')); % Only this newly created private directory.
disk=java.io.File(folder);
assert(double(disk.getUsableSpace)>bytes*1.1,'isolation:DiskSpace','Insufficient space for full waveform export.');
cross=zeros(width); total=zeros(width,1); count=0;
for population=1:2
    if population==1, name='ap'; n=ns; else, name='noise'; n=nn; end
    fid=fopen(fullfile(folder,[name '.bin']),'w','ieee-le');
    assert(fid>=0,'isolation:Export','Cannot create waveform export.');
    closeFile=onCleanup(@()fclose(fid)); timer=tic; last=0;
    for first=1:cfg.FeatureBatchSize:n
        rows=first:min(first+cfg.FeatureBatchSize-1,n);
        if population==1, batch=a(rows,:); else, batch=extractNoise(rows); end
        assert(isequal(size(batch),[numel(rows),width]) && all(isfinite(batch),'all'), ...
            'isolation:Export','Invalid waveform features.');
        assert(fwrite(fid,batch','double')==numel(batch),'isolation:Export','Waveform write failed.');
        cross=cross+batch'*batch; total=total+sum(batch,1)'; count=count+size(batch,1);
        if cfg.Progress && toc(timer)-last>=10
            last=toc(timer); fprintf('  Preparing %s waveforms: %.0f%%\n',name,100*rows(end)/n);
        end
    end
    clear closeFile
end
% PCA supplies orthonormal SEARCH BOUNDS only. Final distances retain all
% 144 original coordinates; no explained-variance cutoff enters the metrics.
covariance=cross/count-(total/count)*(total/count)';
[v,e]=eig((covariance+covariance')/2,'vector'); [~,order]=sort(e,'descend');
q=v(:,order(1:32));
assert(norm(q'*q-eye(32),2)<1e-10,'isolation:Projection','Invalid orthonormal search basis.');
fid=fopen(fullfile(folder,'projection32.bin'),'w','ieee-le'); assert(fid>=0);
closer=onCleanup(@()fclose(fid)); fwrite(fid,q','double'); clear closer
run_native(fullfile(nativeDir,'isolation_knn.exe'), ...
    {folder,folder,num2str(cfg.Threads)},cfg.Progress);
k=jsondecode(fileread(fullfile(folder,'full_scores.json')));
assert(k.FNCount>=0 && k.FNCount<=nn && k.FPCount>=0 && k.FPCount<=ns && ...
    abs(k.FN-k.FNCount/(k.FNCount+ns))<1e-14 && abs(k.FP-k.FPCount/ns)<1e-14, ...
    'isolation:NativeOutput','Native FN/FP counts or normalization failed validation.');
if zeroSpread
    % IS has no positive distance scale, but paper FN/FP remain defined.
    p=nan(ns,1);
    s=struct('APCount',ns,'NoiseCount',nn,'IS',NaN,'CharacteristicDistance',0,'Seconds',0);
else
oldPath=getenv('PATH'); restorePath=onCleanup(@()setenv('PATH',oldPath));
setenv('PATH',[fullfile(cuda,'bin') pathsep oldPath]);
run_native(fullfile(nativeDir,'isolation_is.exe'), ...
    {folder,folder,cfg.ISPrecision,sprintf('%.17g',cfg.DistanceScale)},cfg.Progress);
s=jsondecode(fileread(fullfile(folder,'is_result.json')));
fid=fopen(fullfile(folder,'affinities.bin'),'r','ieee-le'); assert(fid>=0);
closer=onCleanup(@()fclose(fid)); p=fread(fid,Inf,'double=>double'); clear closer
end
assert(k.APCount==ns && k.NoiseCount==nn && s.APCount==ns && s.NoiseCount==nn && ...
    numel(p)==ns && (zeroSpread || (all(isfinite(p)) && all(p>=0 & p<=1) && abs(mean(p)-s.IS)<1e-12)), ...
    'isolation:NativeOutput','Native result counts/affinities failed validation.');
fn=k.FN; fp=k.FP; is=s.IS;
precision='FP32 affinity arithmetic; FP64 scale and sums';
if strcmp(cfg.ISPrecision,'double'), precision='FP64'; end
metric=struct('status','ok','spikes',ns,'noise',nn,'characteristicDistance',s.CharacteristicDistance, ...
    'K',cfg.K,'Completed',true(1,3),'TargetAffinity',p,'FNCount',k.FNCount, ...
    'FPCount',k.FPCount,'FNDefinition','Nfn/(Nfn+Ns)','MetricVersion','Joshua2007-full-v1', ...
    'Backend','native','ISPrecision',precision,'KNNSeconds',k.TotalSeconds,'ISSeconds',s.Seconds);
if zeroSpread, metric.status='zero_spike_spread'; metric.Completed(1)=false; end
end

function run_native(executable,args,progress)
% Direct launch avoids shell interpolation. Blocking line reads do not poll.
% Child processes emit progress about every 10 seconds and at stage ends.
command=java.util.ArrayList; command.add(java.lang.String(executable));
for k=1:numel(args), command.add(java.lang.String(args{k})); end
builder=java.lang.ProcessBuilder(command); builder.redirectErrorStream(true);
environment=builder.environment; environment.put('PATH',getenv('PATH'));
process=builder.start; stop=onCleanup(@()process.destroy);
reader=java.io.BufferedReader(java.io.InputStreamReader(process.getInputStream));
lastLines={};
while true
    line=reader.readLine; if isequal(line,[]), break; end
    text=char(line); if progress, fprintf('%s\n',text); end
    lastLines{end+1}=text; if numel(lastLines)>10, lastLines(1)=[]; end %#ok<AGROW>
end
code=process.waitFor;
assert(code==0,'isolation:NativeFailed','Native calculation failed (%d): %s',code,strjoin(lastLines,newline));
end

function rows=display_rows(n,limit)
rows=unique(round(linspace(1,n,min(n,limit))));
end
