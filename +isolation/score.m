function [result,diagnostics] = score(recording,options)
%SCORE Compute IS, FN, FP and two SNRs for one unit in stable recording periods.
%   [result,diagnostics] = isolation.score(recording,options)
%   Required fields: Signal, SpikeSamples, SampleRateHz, StableIntervals.
%   Optional fields: FileBounds (acquisition joins), ID (label).
%   All intervals/indices are inclusive, one-based samples in Signal.
%   Set options.CaptureQC=true for isolation.plotUnit(diagnostics).
%   Raw data must be an unfiltered single channel. Both AP/noise waveforms
%   receive the same preprocessing described in docs/methods.md.
if nargin<2, options=struct; end
validateattributes(options,{'struct'},{'scalar'});
[rec,mask,id]=normalize_recording(recording);
% A conflicting sampling rate is usually an adapter error; do not hide it.
if isfield(options,'SampleRateHz')
    assert(isequal(double(options.SampleRateHz),double(rec.SampleRateHz)), ...
        'isolation:SampleRateConflict','Use the actual recording.SampleRateHz; options must agree.');
end
[result,~,~,scoring]=score_stable_recording(rec,mask,options);
diagnostics=struct('ID',id,'Result',result,'Scoring',scoring);
end
