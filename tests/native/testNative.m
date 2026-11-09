function tests=testNative
% Optional integration checks; build native executables and set CUDA_PATH.
tests=functiontests(localfunctions);
end

function testNativeAgainstMatlab(t)
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
old=path; restore=onCleanup(@()path(old)); addpath(root);
assumeTrue(t,ispc && isfile(fullfile(root,'native','isolation_knn.exe')) && ...
    isfile(fullfile(root,'native','isolation_is.exe')) && ~isempty(getenv('CUDA_PATH')));
stream=RandStream('mt19937ar','Seed',38); a=randn(stream,38,144); b=randn(stream,67,144)+.2;
[ref,rd]=isolation.scoreWaveforms(a,b,struct('Backend','matlab','Progress',false));
for precision={'mixed','double'}
    options=struct('Backend','native','Progress',false,'ISPrecision',precision{1});
    [r,d]=isolation.scoreWaveforms(a,b,options);
    verifyEqual(t,[r.fnScore r.fpScore],[ref.fnScore ref.fpScore]);
    tolerance=1e-6; if strcmp(precision{1},'double'), tolerance=1e-12; end
    verifyEqual(t,d.Scoring.Metric.TargetAffinity,rd.Scoring.Metric.TargetAffinity,'AbsTol',tolerance);
end
end
