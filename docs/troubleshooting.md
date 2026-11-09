# Troubleshooting

Always inspect `result.Status`, `diagnostics.Scoring.Warnings`, population
counts and `MetricStatistics`. A batch exception instead has `ErrorIdentifier`
and `ErrorMessage` directly in that unit's diagnostics.

| Status/error | Meaning and next check |
| --- | --- |
| `ok` | Metrics calculated; inspect biological QC separately |
| `ok_low_population` | A full population is below the reporting threshold; no downsampling occurred |
| `no_stable_trials` | No retained interval was supplied |
| `insufficient_eligible_spikes` | Fewer than two targets survived interval/edge eligibility |
| `flat_template` | No nonzero target phase could be found |
| `invalid_noise_threshold` | Weak-tail threshold is nonnegative/nonfinite; inspect polarity/alignment and spike labels |
| `insufficient_events` | Fewer than two APs or no detected noise; IS/FN/FP are undefined here |
| `insufficient_knn_population` | Fewer than K+1 events; FN/FP are NaN while IS may be valid |
| `zero_spike_spread` | Identical target features give no positive IS scale; FN/FP may remain valid |
| `waveform_memory_limit` | All-AP SNR finished, but AP cache exceeds the workspace guard |
| `isolation:MatlabBudget` | Reference pair budget exceeded; build native or deliberately raise the guard |
| `isolation:NativeShape` | Native platform/feature-width/K constraints are not met |
| `isolation:NativeMissing` | Build the optional executables in `native/` |
| `isolation:CUDA` / `isolation:NativeFailed` | Check toolkit/runtime, driver, GPU availability, memory and the included native error text |
| `isolation:MemoryBudget` / `isolation:DiskSpace` | Required full populations exceed configured/available resources |
| `isolation:SpikeIndices` | Duplicate, unsorted or out-of-range spike indices; fix the adapter with provenance |

Raw results additionally report `InputSpikes`, `EligibleSpikes`, `EligibleNoise`,
`InputSamples`, `NumAPsSNR`, and `NumBaselineSamples`. `NumAPs`/`NumNoise` count
events entering metric scoring; they can remain zero when preprocessing fails
earlier. Use the eligibility counts to identify where events were lost.

An empty or NaN metric is not zero. Do not lower the noise threshold merely to
make a status disappear: doing so changes the comparison population.

For apparently missing amplitude points, green spans show retained stable
intervals and gray spans show excluded time. A green span can still contain no
eligible spikes. Dotted acquisition boundaries do not imply instability.

For slow jobs, start with one representative unit, inspect populations, and
choose the appropriate backend. Plot subsets and `DisplayEvents` do not limit
metric populations. The native backend has explicit prerequisites; automatic
selection does not mean automatic installation or unlimited memory.
