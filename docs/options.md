# Options

`isolation.defaults` returns all current calculation settings. Pass only the
overrides you need; unknown names are rejected. For raw scoring the actual
rate is `recording.SampleRateHz`; a supplied option with a different rate is an
error. Pass settings directly to the scoring function in an options struct.

## Scientific settings

| Option | Default | Purpose |
| --- | --- | --- |
| `BinSamples` | 36 | Cropped waveform length before interpolation |
| `SplitRate` | 4 | Interpolation factor; feature width is their product |
| `K` | 31 | Odd pooled-neighbor count for FN/FP |
| `DistanceScale` | 0.1 | Scale multiplier for the IS characteristic distance |
| `ThresholdFraction` | 0.02 | Weakest fraction of all APs for threshold estimation |
| `ThresholdRatio` | 2 | Weak-tail mean is divided by this value |
| `Polarity` | `'auto'` | `'auto'`, `'negative'` or `'positive'` |
| `ShiftSearchMs` | 3.6 | Shared template-offset search radius |
| `AlignmentSearchMs` | 4/30 | Local integer-trough search radius in ms |
| `FilterPaddingMs` | 50 | Same-file context on each side of a stable segment |
| `FileEdgeGuardMs` | 20 | Guard at acquisition file edges |
| `MinimumReportedSpikes` | 100 | Low-population flag for either AP or noise |
| `SampleRateHz` | 30000 | Required to agree with raw recording when specified |

Change scientific settings deliberately and record them with your results.
Native supports feature width 144 and K=31. Changing the detection threshold
changes which events enter the noise population.

## Execution and diagnostics

| Option | Default | Purpose |
| --- | --- | --- |
| `Backend` | `'auto'` | `'matlab'`, `'native'`, or selection by pair budget |
| `ISPrecision` | `'mixed'` | Native IS only: `'mixed'` or `'double'` |
| `Threads` | 12 | Native CPU worker threads; lower on smaller/shared machines |
| `MaxMemoryMB` | 32768 | Scoring workspace guard; not a total process RAM limit |
| `MaxMatlabPairs` | 5e7 | Guard on `Ns²+Ns*Nn+Nn²`; also controls `auto` |
| `TemporaryDirectory` | `tempdir` | Existing folder for native waveform exports |
| `CUDARoot` | `getenv('CUDA_PATH')` | CUDA toolkit root for native runtime libraries |
| `DistanceBlockSize` | 256 | MATLAB query block size |
| `DistanceColumnBlockSize` | 2048 | MATLAB reference block size |
| `FeatureBatchSize` | 4096 | Preprocessing/export batch size |
| `CaptureQC` | false | Collect diagnostics for unit plotting (raw API only) |
| `DisplayEvents` | 3000 | Maximum displayed/stored waves per population |
| `Progress` | true | Stage progress; native updates roughly every 10 seconds |

`auto` uses MATLAB within the pair budget and attempts native above it. On
systems without native support, explicitly select MATLAB and increase `MaxMatlabPairs`
only after assessing runtime. `MaxMemoryMB` excludes raw recording storage and
is a guard, not an operating-system memory reservation.

Native exports require approximately `8*(Ns+Nn)*D` bytes of free temporary disk
plus headroom. Estimated native scoring RAM is guarded at three times those
waveform bytes. A private temporary subdirectory is removed on normal success
or MATLAB exceptions; forcibly killing MATLAB can leave it behind.

For a smaller machine, reduce `Threads`/`MaxMemoryMB`, choose a disk with enough
space, and first validate one representative unit. Reducing `DisplayEvents`
does **not** reduce scoring cost.

## Batch-only options

| Option | Default | Purpose |
| --- | --- | --- |
| `Score` | `struct` | Calculation overrides from tables above |
| `Progress` | true | Unit count, rough ETA, elapsed time and status |
| `Recompute` | false | Recalculate complete saved results when true |
| `ErrorMode` | `'continue'` | Retain errors and continue, or `'error'` to stop |

Set `jobs(k).CaptureQC` per unit instead of `options.Score.CaptureQC`. Batch
progress is based on completed jobs. The first
unit has no ETA. Job duration includes loading/preprocessing. Populations can
vary markedly between units, so the estimate may change substantially.
