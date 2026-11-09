function options = defaults(overrides)
%DEFAULTS Return calculation settings, optionally applying validated overrides.
%   options = isolation.defaults(struct('Backend','matlab','Progress',false));
%   See docs/options.md. DisplayEvents controls plotted waveform subsets.
%   SampleRateHz comes from recording in isolation.score.
if nargin<1, overrides=struct; end
validateattributes(overrides,{'struct'},{'scalar'});
options=isolation_config(overrides);
end
