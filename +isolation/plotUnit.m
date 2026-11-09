function fig = plotUnit(diagnostics,visible)
%PLOTUNIT Display saved unit QC without loading data or recalculating scores.
%   Set CaptureQC=true when calling isolation.score, or on selected batch
%   jobs. Missing QC produces an explanatory panel. visible='off' is useful
%   with exportgraphics(fig,'unit.png','Resolution',150).
if nargin<2, visible='on'; end
validateattributes(diagnostics,{'struct'},{'scalar'});
assert(all(isfield(diagnostics,{'ID','Result'})),'isolation:Diagnostics','Expected diagnostics returned by score or batch.');
visible=validatestring(visible,{'on','off'});
fig=plot_isolation_unit(diagnostics,visible);
end
