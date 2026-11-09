function tests=testAPI
% Public contract tests use synthetic data only; no recording files needed.
tests=functiontests(localfunctions);
end

function setupOnce(t)
t.TestData.OldPath=path;
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); addpath(root);
t.TestData.Record=isolation.demoData(7);
t.TestData.Options=struct('Backend','matlab','Progress',false,'CaptureQC',true);
[t.TestData.Result,t.TestData.Diagnostics]=isolation.score(t.TestData.Record,t.TestData.Options);
end

function teardownOnce(t)
path(t.TestData.OldPath);
end

function testStableGapsAndFileJoins(t)
r=t.TestData.Record; d=t.TestData.Diagnostics.Scoring;
verifyEqual(t,d.StableBlocks,r.StableIntervals);
verifyEqual(t,d.QC.CrossBoundaryWindows,0);
verifyTrue(t,d.QC.CountsMatch);
times=d.QC.TargetTimeSeconds*r.SampleRateHz+1;
verifyTrue(t,all(times<=r.StableIntervals(1,2) | times>=r.StableIntervals(2,1)));
r.StableIntervals=[1001 179000];
[~,joined]=isolation.score(r,t.TestData.Options);
verifyEqual(t,joined.Scoring.Segments,[1001 90000 1;90001 179000 2]);
end

function testDisplayCapDoesNotSampleMetrics(t)
options=t.TestData.Options; options.DisplayEvents=12; options.CaptureQC=false;
r=isolation.score(t.TestData.Record,options);
verifyEqual(t,r,t.TestData.Result);
end

function testSignInvariance(t)
record=t.TestData.Record; record.Signal=-record.Signal;
[r,d]=isolation.score(record,t.TestData.Options);
verifyEqual(t,r,t.TestData.Result);
verifyEqual(t,d.Scoring.PolarityMultiplier,-t.TestData.Diagnostics.Scoring.PolarityMultiplier);
end

function testRNGPreserved(t)
state=rng; restore=onCleanup(@()rng(state));
isolation.demoData(12); verifyEqual(t,rng,state);
isolation.score(t.TestData.Record,struct('Backend','matlab','Progress',false));
verifyEqual(t,rng,state);
record=t.TestData.Record; record.SpikeSamples=record.SpikeSamples(:);
[result,d]=isolation.score(record,t.TestData.Options);
verifyEqual(t,result,t.TestData.Result);
verifyEqual(t,d,t.TestData.Diagnostics); % Row/column inputs canonicalize alike.
end

function testRejectMalformedInputs(t)
r=t.TestData.Record; r.SpikeSamples=[100 100];
verifyError(t,@()isolation.score(r),'isolation:SpikeIndices');
r.SpikeSamples=[200 100]; verifyError(t,@()isolation.score(r),'isolation:SpikeIndices');
r=t.TestData.Record; r.StableIntervals=[1 20;15 30];
verifyError(t,@()isolation.score(r),'isolation:Intervals');
r=t.TestData.Record; r.FileBounds=[1 10;12 numel(r.Signal)];
verifyError(t,@()isolation.score(r),'isolation:FileCoverage');
r=t.TestData.Record; r.SpikeTimes=r.SpikeSamples;
verifyError(t,@()isolation.score(r),'isolation:RecordingFields');
verifyError(t,@()isolation.score(t.TestData.Record,struct('SampleRateHz',20000)), ...
    'isolation:SampleRateConflict');
verifyError(t,@()isolation.defaults(struct('NumWaveforms',5000)),'isolation:UnknownOption');
end

function testEmptyAndShortPopulations(t)
r=t.TestData.Record; r.StableIntervals=zeros(0,2);
result=isolation.score(r); verifyEqual(t,result.Status,'no_stable_trials');
r=t.TestData.Record; r.SpikeSamples=[];
result=isolation.score(r); verifyEqual(t,result.Status,'insufficient_eligible_spikes');
r=t.TestData.Record; r.StableIntervals=[5000 5001];
result=isolation.score(r); verifyEqual(t,result.Status,'insufficient_eligible_spikes');
end

function testAdjacentIntervalsMerge(t)
r=t.TestData.Record; r.StableIntervals=[1001 40000;40001 72000;108001 179000];
[result,d]=isolation.score(r,t.TestData.Options);
verifyEqual(t,d.Scoring.StableBlocks,t.TestData.Record.StableIntervals);
verifyEqual(t,result,t.TestData.Result);
end

function testNonfiniteSelectedData(t)
r=t.TestData.Record; r.Signal(5000)=NaN;
verifyError(t,@()isolation.score(r),'isolation:NonfiniteSignal');
% Missing values wholly outside selected intervals AND filtering padding are
% harmless; this is deliberately different from a NaN near a selected block.
r=t.TestData.Record; r.Signal(90000)=NaN;
result=isolation.score(r,struct('Backend','matlab','Progress',false));
verifyEqual(t,result,t.TestData.Result);
end

function testResumeDoesNotCallLoader(t)
jobs=struct('ID',"saved",'Load',@mustNotLoad,'Result',t.TestData.Result);
[r,d]=isolation.batch(jobs,struct('Progress',false));
verifyEqual(t,r{1},t.TestData.Result); verifyTrue(t,d{1}.ReusedExisting);
jobs.Result.IS=NaN;
verifyError(t,@()isolation.batch(jobs,struct('Progress',false,'ErrorMode','error')), ...
    'test:LoaderWasCalled');
jobs.Result=t.TestData.Result;
verifyError(t,@()isolation.batch(jobs,struct('Progress',false,'ErrorMode','error','Recompute',true)), ...
    'test:LoaderWasCalled');
end

function testBatchMatchesSingleAndRetainsErrors(t)
jobs(1)=struct('ID',"session-a/unit-1",'Load',@()t.TestData.Record,'CaptureQC',true);
jobs(2)=struct('ID',"session-b/unit-1",'Load',@mustNotLoad,'CaptureQC',false);
old=warning('off','isolation:JobFailed'); restore=onCleanup(@()warning(old));
[r,d]=isolation.batch(jobs,struct('Progress',false,'Score',struct('Backend','matlab')));
verifyEqual(t,r{1},t.TestData.Result); verifyEqual(t,r{2}.Status,'error');
verifyEqual(t,d{2}.ErrorIdentifier,'test:LoaderWasCalled');
verifyEqual(t,d{1}.ID,"session-a/unit-1");
verifyTrue(t,isfield(d{1}.Scoring,'QC'));
end

function testPlottingDoesNotChangeResults(t)
before=t.TestData.Diagnostics;
fig=isolation.plotUnit(before,'off'); c=onCleanup(@()close(fig));
verifyEqual(t,numel(findall(fig,'Type','axes')),10); % Nine panels + PCA inset.
verifyEqual(t,before,t.TestData.Diagnostics);
fig2=isolation.plotSummary({t.TestData.Result},'Synthetic test','off');
c2=onCleanup(@()close(fig2));
verifyEqual(t,numel(findall(fig2,'Type','axes')),9);
end

function record=mustNotLoad
record=[]; %#ok<NASGU> Deliberately fail after defining the output for analysis.
error('test:LoaderWasCalled','This deliberate failing loader must be skipped on resume.');
end
