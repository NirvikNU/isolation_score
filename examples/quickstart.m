% Run from any folder; only the package root belongs on the MATLAB path.
packageRoot=fileparts(fileparts(mfilename('fullpath')));
addpath(packageRoot);
recording=isolation.demoData;
options=struct('Backend','matlab','CaptureQC',true);
[result,diagnostics]=isolation.score(recording,options);
disp(result);
isolation.plotUnit(diagnostics);
isolation.plotSummary({result},'Synthetic demonstration');
% Optional: save('my_results.mat','result','diagnostics');
