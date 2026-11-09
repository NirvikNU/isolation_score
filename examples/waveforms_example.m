% Prepared waveform API: both populations must share exactly the same
% preprocessing, units, alignment and coordinate grid. No raw filtering is
% performed here. Real noise events must come from the chosen detector.
addpath(fileparts(fileparts(mfilename('fullpath'))));
stream=RandStream('mt19937ar','Seed',4);
grid=linspace(-1,2,144); template=-exp(-grid.^2/.08);
ap=template+.15*randn(stream,100,144);
noise=.45*template+.25*randn(stream,200,144);
ap=ap-mean(ap,2); noise=noise-mean(noise,2);
[result,diagnostics]=isolation.scoreWaveforms(ap,noise,struct('Backend','matlab'));
disp(result); % snrBeforeAp is NaN because no raw baseline was supplied.
