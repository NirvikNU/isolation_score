function recording = demoData(seed)
%DEMODATA Small synthetic one-channel recording; does not alter caller RNG.
%   Background and overlapping target/non-target waveforms are illustrative,
%   not a biological validation dataset. Stable periods contain a real gap.
if nargin<1, seed=7; end
validateattributes(seed,{'numeric'},{'real','finite','scalar','integer','nonnegative'});
stream=RandStream('mt19937ar','Seed',seed);
fs=30000; n=180000; raw=.045*randn(stream,1,n); q=-30:30;
spikes=3000:1400:n-3000;
spikes=sort([spikes,spikes(1:6:end)+150]); % A few 5-ms pairs make ISI QC visible.
target=-exp(-q.^2/18)+.25*exp(-(q-10).^2/25);
for k=1:numel(spikes)
    raw(spikes(k)+q)=raw(spikes(k)+q)+(1+.1*randn(stream))*target;
end
other=3500:1900:n-3000;
for k=1:numel(other)
    raw(other(k)+q)=raw(other(k)+q)+(.7+.12*randn(stream))* ...
        (-exp(-q.^2/30)+.3*exp(-(q-12).^2/30));
end
recording=struct('Signal',raw,'SpikeSamples',spikes,'SampleRateHz',fs, ...
    'StableIntervals',[1001 72000;108001 179000], ...
    'FileBounds',[1 90000;90001 n],'ID',"synthetic unit");
end
