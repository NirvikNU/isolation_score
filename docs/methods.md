# Calculation definitions and implementation choices

The methodological reference is [Joshua et al. (2007)](https://doi.org/10.1016/j.jneumeth.2007.03.012).
The description below defines this implementation, including choices that
should be reported alongside scores. This is not a claim that every original
acquisition or preprocessing detail is reproduced.

## Populations and fixed preprocessing

1. Merge physically adjacent stable intervals, retain gaps, and split at every
   acquisition-file boundary. High-pass filter each segment with same-file
   padding: Butterworth design from a 300-Hz passband, 100-Hz stopband, 3-dB
   passband loss and 10-dB stopband attenuation, applied using `filtfilt`.
2. Reject edge events that cannot support the complete waveform and alignment
   searches. Default acquisition-edge guard is 20 ms; filter padding is 50 ms.
3. Use **all eligible APs** to form the initial mean template. Choose its
   dominant positive or negative phase once per unit, normalize the entire
   filtered signal to a negative target phase, and determine the shared global
   timing offset. `Polarity` can explicitly select positive/negative instead.
4. Apply local integer trough alignment and identical spline interpolation,
   sub-sample trough alignment, crop and per-waveform mean subtraction to AP
   and noise. Default 36 samples × interpolation factor 4 gives 144 columns;
   at 30 kHz the cropped waveform spans 1.2 ms.
5. At the all-AP mean waveform's trough column, sort AP values. Average the
   weakest `ceil(0.02*Ns)` (least-negative) values and divide by 2. This
   negative value is the noise-detection threshold. The weakest count is at
   least one. This **2% and 1/2 convention is an explicit implementation
   setting**, not an adaptive event-ratio rule.
6. Group samples below threshold into excursions (default grouping gap scales
   from four samples at 30 kHz). Take the first minimum in each group and apply
   the same waveform processing. Exclude target centers within a timing
   tolerance scaling from ±2 samples at 30 kHz. Threshold events from other
   units intentionally remain in the background population.

All eligible events are retained. There is no incremental sampling, uncertainty
stopping rule, or reduction to a fixed AP/noise ratio. Very few noise events
are reported rather than repaired by changing the threshold. A nonnegative
threshold is invalid and must be investigated.

## Isolation score

Let `A` be the prepared target waveforms, `B` the prepared noise waveforms,
`Ns` the number of APs, and `d(x,y)` the full-waveform Euclidean distance.

```text
h = DistanceScale * sum(i=1..Ns, j=1..Ns) d(Ai,Aj) / Ns^2
w(x,y) = exp(-d(x,y) / h)
pi = sum(j != i) w(Ai,Aj)
     / [sum(j != i) w(Ai,Aj) + sum(k=1..Nn) w(Ai,Bk)]
IS = mean(pi)
```

`DistanceScale=0.1` by default. Zero self-distances enter the `Ns²`
normalization of `h`; self-affinity is excluded from `pi`. The kernel uses
Euclidean distance, not squared distance. AP–AP and AP–noise terms are needed
for IS; noise–noise terms are not. Higher IS indicates stronger relative
affinity to the target population under this preprocessing and detector.

## False-negative and false-positive scores

Classify each event by the majority among its **K=31 nearest other events** in
the pooled AP/noise population. Every distance uses all waveform coordinates.
Self is excluded. Exact distance ties give APs priority over noise; within a
population, original order breaks remaining ties. K is never silently reduced.

```text
Nfn = number of noise events classified as AP
Nfp = number of AP events classified as noise
fnScore = Nfn / (Nfn + Ns)
fpScore = Nfp / Ns
```

These use the paper's FN/FP normalizations. They are classifier estimates
relative to the supplied labels and background detector, not measured errors
against independently known spike identities. Noise–noise comparisons are
needed to determine noise-query votes; the implementation avoids unnecessary
comparisons using certified bounds or early rejection.

## Two SNRs

Let `mu` be the mean of **all** prepared AP waveforms and `p2p=max(mu)-min(mu)`.

- `snrAP = p2p / std((A-mu)(:))`. This pools waveform residual elements and uses
  the sample-SD denominator `Ns*D-1`.
- `snrBeforeAp = p2p / std(baselineSamples)`. Baseline windows are one
  `BinSamples` interval ending one bin before eligible isolated target spikes.
  The target must be separated from the preceding target by more than three
  bins; the first target in a segment has no preceding interval and is omitted.
  Baselines cannot cross a stable/file boundary.

There is **no division by five** in either measure. Thus report these SNR
definitions explicitly when comparing with other implementations. Undefined
zero-variance baselines or unavailable pre-spike data yield NaN.

## Computation and precision

The MATLAB reference computes blocked FP64 Euclidean distances and maintains
affinity ratios with a shared exponent offset for numerical stability. It
stores no full all-event distance matrix. A conservative pair-count guard
prevents accidentally launching impractical reference calculations.

The optional native backend computes FN/FP on an AVX2 CPU using 32 orthonormal
PCA coordinates as **search bounds only**. Final candidate distances retain all
144 original coordinates. It never substitutes a reduced-dimensional metric.
GPU IS includes all AP–AP and AP–noise pairs and streams noise chunks. The
characteristic distance and accumulated sums use FP64. Default `mixed` uses
FP32 dot products/affinity arithmetic with cancellation safeguards; `double`
uses FP64 affinity arithmetic. Mixed precision is a numerical approximation,
not statistical sampling. Record backend and precision from `MetricStatistics`.

## Interpretation and limitations

- At least two APs and one noise event are required here. FN/FP additionally
  need at least `K+1` pooled events. Zero AP spread makes IS undefined while
  classifier scores can still be calculated.
- `MinimumReportedSpikes=100` flags either population below 100. This is a QC
  heuristic, **not** a minimum population prescribed by the paper or a proof
  of estimator reliability. Full population scoring does not eliminate
  detector/label bias or finite-recording uncertainty.
- Scores depend on event selection, stable periods, feature grid and noise
  definition. Changing those settings changes the quantity being measured.
- Polarity is fixed per unit, not optimized to maximize IS. Nearly balanced
  biphasic templates receive a diagnostic warning.
- The default fixed sample width was validated at 30 kHz. Different rates,
  filtering or waveform durations require appropriate scientific validation.
- SNR, IS, FN and FP capture different properties. Inspect waveforms, stability,
  ISIs and calculation statuses together rather than imposing a universal
  quality cutoff from this package alone.
