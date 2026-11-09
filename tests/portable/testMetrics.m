function tests=testMetrics
% Compare with an independent, exhaustive implementation on small inputs.
tests=functiontests(localfunctions);
end

function setupOnce(t)
t.TestData.OldPath=path;
addpath(fileparts(fileparts(fileparts(mfilename('fullpath')))));
stream=RandStream('mt19937ar','Seed',94);
t.TestData.AP=randn(stream,38,144); t.TestData.Noise=randn(stream,67,144)+.2;
t.TestData.Options=struct('Backend','matlab','Progress',false);
end

function teardownOnce(t)
path(t.TestData.OldPath);
end

function testFullMetricsAgainstExhaustive(t)
a=t.TestData.AP; b=t.TestData.Noise;
[r,d]=isolation.scoreWaveforms(a,b,t.TestData.Options);
[fn,fp,score,p]=reference(a,b,31,.1);
verifyEqual(t,[r.fnScore r.fpScore],[fn fp]);
verifyEqual(t,r.IS,score,'AbsTol',1e-12);
verifyEqual(t,d.Scoring.Metric.TargetAffinity,p,'AbsTol',1e-12);
average=mean(a,1); residual=a-average;
verifyEqual(t,r.snrAP,(max(average)-min(average))/std(residual(:)),'AbsTol',1e-12);
verifyTrue(t,isnan(r.snrBeforeAp));
end

function testRareNoiseAndDistanceTies(t)
a=t.TestData.AP;
for n=[1 2 15 16 17 31]
    b=repmat(a(1,:),n,1);
    [fn,fp]=reference(a,b,31,.1);
    r=isolation.scoreWaveforms(a,b,t.TestData.Options);
    verifyEqual(t,[r.fnScore r.fpScore],[fn fp]);
    verifyTrue(t,isfinite(r.IS)); verifyEqual(t,r.NumNoise,n);
end
end

function testDegeneratePopulations(t)
options=t.TestData.Options; a=t.TestData.AP; b=t.TestData.Noise;
r=isolation.scoreWaveforms(a,zeros(0,144),options);
verifyTrue(t,isnan(r.IS)); verifyEqual(t,r.Status,'insufficient_events');
r=isolation.scoreWaveforms(ones(38,144),b,options);
verifyTrue(t,isnan(r.IS)); verifyTrue(t,isfinite(r.fnScore) && isfinite(r.fpScore));
verifyEqual(t,r.Status,'zero_spike_spread');
r=isolation.scoreWaveforms(a(1:2,:),b(1:2,:),options);
verifyTrue(t,isfinite(r.IS)); verifyTrue(t,isnan(r.fnScore) && isnan(r.fpScore));
verifyEqual(t,r.Status,'insufficient_knn_population');
end

function testResourceGuardAndValidation(t)
options=t.TestData.Options; options.MaxMatlabPairs=1;
verifyError(t,@()isolation.scoreWaveforms(t.TestData.AP,t.TestData.Noise,options), ...
    'isolation:MatlabBudget');
verifyError(t,@()isolation.scoreWaveforms(t.TestData.AP,ones(20,12)), ...
    'isolation:FeatureWidth');
options=t.TestData.Options; options.CaptureQC=true;
verifyError(t,@()isolation.scoreWaveforms(t.TestData.AP,t.TestData.Noise,options), ...
    'isolation:WaveformQC');
end

function [fn,fp,score,p]=reference(a,b,k,scale)
x=[a;b]; ns=size(a,1); n=size(x,1); distances=zeros(n);
for i=1:n, distances(i,:)=sqrt(sum((x-x(i,:)).^2,2)); end
h=scale*sum(distances(1:ns,1:ns),'all')/ns^2;
p=zeros(ns,1);
for i=1:ns
    weights=exp(-distances(i,:)/h); weights(i)=0;
    p(i)=sum(weights(1:ns))/sum(weights);
end
score=mean(p); labels=false(n,1);
for i=1:n
    row=distances(i,:)'; row(i)=Inf;
    [~,order]=sortrows([row,(1:n)'],[1 2]);
    labels(i)=sum(order(1:k)<=ns)>k/2;
end
missed=sum(labels(ns+1:end)); fn=missed/(missed+ns); fp=mean(~labels(1:ns));
end
