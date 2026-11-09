# Input contract

## Raw recording: `isolation.score`

One call represents one target unit on one channel. Several units can share the
same channel. The adapter is responsible for selecting the correct channel and
converting timestamps to the signal's sample coordinate system.

| Field | Required | Format |
| --- | --- | --- |
| `Signal` | yes | Real, nonsparse numeric vector; one unfiltered channel |
| `SampleRateHz` | yes | Actual sampling rate, finite and >600 Hz |
| `SpikeSamples` | yes | Sorted, unique, positive integer target-spike indices |
| `StableIntervals` | yes | N×2 sorted, disjoint inclusive `[first last]` sample intervals |
| `FileBounds` | no | N×2 inclusive acquisition boundaries; default `[1 numel(Signal)]` |
| `ID` | no | Text or numeric scalar label; default `"unit"` |

Unknown fields/options are rejected to catch spelling mistakes. Single/double
and integer raw signal arrays are accepted; filtering operates in double
precision. Amplitude units must be consistent within the recording; they can
be microvolts or ADC units. Do not scale AP and noise differently.

All indices are **one-based**. A timestamp `t` seconds after the first sample
maps to `round(t*SampleRateHz)+1` only if the acquisition's time origin and
rounding convention support that conversion. Prefer original integer sample
indices. For concatenated files, add the number of samples in preceding files
to each file-local index. Validate units and origin before conversion.

The input signal retains unstable samples. `StableIntervals` selects periods
for scoring without changing their original coordinates. Adjacent intervals
merge; separated intervals remain separate. `zeros(0,2)` means no stable
periods, producing `Status='no_stable_trials'`; it does not select everything.
To select a whole recording, explicitly use `[1 numel(Signal)]`.

For trials 1–100 and 150–400 marked stable, construct two intervals:

```matlab
stable = [trialStart(1) trialEnd(100); ...
          trialStart(150) trialEnd(400)];
```

Do not enter every trial separately if the intervals between those trials
should be retained. Trials or behavioral events are not otherwise required.

`FileBounds` must partition the entire signal without uncovered samples or
overlap. Every acquisition restart or discontinuous concatenation must be
represented, even within a stable interval. Time axes use concatenated sample
time; they cannot reconstruct wall-clock gaps between files.

The filter uses real padding from the same acquisition file, including context
outside stable periods. Event centers and complete waveform/search windows
must remain in an eligible stable/file segment. Thus some spikes near edges
are intentionally excluded. Selected signal/filter context must be finite;
NaNs far outside all selected intervals and padding are never processed.

Do not concatenate extracted stable pieces into an apparently continuous
signal. Do not run filtering across a file join. The package needs the
boundaries to enforce those rules itself.

## Prepared waveforms: `isolation.scoreWaveforms`

`ap` is `Ns × D`; `noise` is `Nn × D`. Every row is one event. Both matrices
must be real, finite, nonsparse and use the same coordinates and units.
An empty population must retain its width, for example `zeros(0,144)`.

This lower-level entry point deliberately performs no filtering, interpolation,
alignment, polarity choice, per-event centering or threshold detection. The
caller must apply the same preprocessing to both populations. For equivalence
with the raw-data pipeline, follow [Methods](methods.md) exactly, including
event selection. Do not duplicate or oversample waveforms to change the class
ratio: the ratio is part of the metric's input population.

MATLAB supports arbitrary positive feature width; native acceleration requires
144 columns. The raw pipeline produces `BinSamples × SplitRate` columns (144
by default). A different sample rate changes the physical duration of a fixed
`BinSamples`; choose and report settings deliberately.

Pre-spike SNR is NaN because raw baseline windows are absent. Waveform SNR uses
all supplied APs. Diagnostic raw traces/stable timelines are unavailable for
this API; `CaptureQC=true` is rejected rather than suggesting those checks ran.

## Batch jobs and saved results

Each job has a distinct `ID` and either `Recording` or `Load`. A loader is a
zero-argument function handle returning the recording struct. Optional fields
are `CaptureQC`, `Result`, and `Diagnostics`. Outputs are ordered cell arrays.

Units are processed sequentially; each loader runs at most once per attempted
job. The package does not discover sessions or cache acquisition channels on
your behalf. Reuse a loaded signal in your adapter when units share it, or use
lazy loaders to limit memory. MATLAB's copy-on-write can share an unchanged
signal, but modifying a shared array may allocate a copy.

Complete finite `IS`, `snrAP`, `snrBeforeAp`, `fnScore` and `fpScore` values
qualify for resume. Zero is valid; NaN/Inf/empty is incomplete. Results are
preserved rather than silently upgraded to a different calculation version.
Resume cannot infer that a changed signal, stable interval or setting should
invalidate an old result; use `Recompute=true` in that case.

`ErrorMode='continue'` records an error result plus identifier/message and
continues to the next job. `'error'` raises immediately. Save outputs explicitly;
batch does not overwrite files or implement automatic disk checkpoints.
