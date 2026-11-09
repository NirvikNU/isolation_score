function [result,diagnostics] = scoreWaveforms(ap,noise,options)
%SCOREWAVEFORMS Compute full-population metrics for already-prepared waveforms.
%   Rows are events; columns are identical waveform coordinates in both
%   populations. This entry point does NOT filter, align, normalize or select
%   noise events. See docs/data-format.md before using externally cut waves.
%   SNRbeforeAp is unavailable (NaN): pre-spike raw samples were not supplied.
%   CaptureQC is not supported here; use isolation.score for raw-data plots.
if nargin<3, options=struct; end
cfg=isolation.defaults(options);
validateattributes(ap,{'numeric'},{'real','finite','2d','nonsparse'});
validateattributes(noise,{'numeric'},{'real','finite','2d','nonsparse'});
assert(size(ap,2)>0 && size(ap,2)==size(noise,2),'isolation:FeatureWidth', ...
    'AP and noise must have the same positive number of waveform columns.');
assert(~cfg.CaptureQC,'isolation:WaveformQC','Use isolation.score to capture raw-data QC.');
ap=double(ap); noise=double(noise);
[statistics,metric]=isolation_full_scores(ap,size(noise,1),@(rows)noise(rows,:),cfg);
snr=NaN;
if size(ap,1)>=2
    average=mean(ap,1); residual=ap-average; sd=std(residual(:));
    if sd>0, snr=(max(average)-min(average))/sd; end
end
status=metric.status;
if strcmp(status,'ok') && min(size(ap,1),size(noise,1))<cfg.MinimumReportedSpikes
    status='ok_low_population';
end
result=struct('IS',statistics.IS.Mean,'snrAP',snr,'snrBeforeAp',NaN, ...
    'fnScore',statistics.fnScore.Mean,'fpScore',statistics.fpScore.Mean, ...
    'NumAPs',size(ap,1),'NumNoise',size(noise,1),'NumAPsSNR',size(ap,1), ...
    'Status',status,'MetricVersion','Joshua2007-full-v1','MetricStatistics',statistics);
diagnostics=struct('ID',"prepared waveforms",'Result',result, ...
    'Scoring',struct('Config',cfg,'Metric',metric,'Warnings',{{ ...
    'Noise population and preprocessing were supplied by the caller; pre-spike SNR is unavailable.'}}));
end
