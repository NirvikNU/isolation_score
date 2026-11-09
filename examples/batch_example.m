% Replace demoData loaders with your own zero-argument recording loaders.
% A loader is called only when its result needs calculating. The job ID can
% identify any session/unit combination; the package assumes no file format.
addpath(fileparts(fileparts(mfilename('fullpath'))));
jobs(1)=struct('ID',"session-A/unit-1",'Load',@()isolation.demoData(7),'CaptureQC',true);
jobs(2)=struct('ID',"session-B/unit-2",'Load',@()isolation.demoData(9),'CaptureQC',false);
options=struct('Score',struct('Backend','matlab'),'Progress',true);
[results,diagnostics]=isolation.batch(jobs,options);
isolation.plotUnit(diagnostics{1});              % Selected unit only.
isolation.plotSummary(results,'All sessions');   % Or results([1 3 5]), etc.
% Resume: attach saved results/diagnostics before calling batch again.
for k=1:numel(jobs)
    jobs(k).Result=results{k}; jobs(k).Diagnostics=diagnostics{k};
end
[resumed,reused]=isolation.batch(jobs,options);
