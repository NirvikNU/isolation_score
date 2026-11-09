# Validation

## Tested configurations

The package is tested on Windows using MATLAB R2024a and Signal Processing
Toolbox. The optional backend was built with Visual Studio
2022, Windows SDK 10.0.22000.0 and CUDA 11.2, and tested on an AVX2 CPU with an
NVIDIA RTX 3080. These are tested configurations, not universal compatibility
claims.

The portable suite contains 15 tests covering:

- Separate stable runs, acquisition splits, adjacent intervals, empty periods,
  short segments, and nonfinite data inside/outside selected filter context.
- Rejection of duplicate/unsorted spikes, overlapping intervals, uncovered file
  coordinates, unknown fields and conflicting sampling rates.
- Exact invariance of all metrics to the display cap and sign inversion.
- Preservation of the caller's random-number state.
- Full metrics against an independent exhaustive Euclidean-distance and
  nearest-neighbor implementation, including distance ties and rare noise.
- No-noise, zero-spread and insufficient-KNN populations, plus resource guards.
- Lazy-loader resume, deliberate recomputation, batch errors and plot creation.

The optional native integration test checks exact agreement of FN/FP with the
MATLAB reference and per-AP IS affinities within `1e-6` for mixed precision and
`1e-12` for double precision on its synthetic fixture. These tolerances are
test criteria on that fixture, not universal error bounds. Native tests skip
when the required executable/platform/toolkit environment is absent.

All three examples use generated data. The README figure is generated from
`isolation.demoData`. This small example runs quickly and illustrates the
low-population status as well as the diagnostic plots.

## Validating a new dataset

Synthetic regressions do not validate a new acquisition adapter's timestamp
origin, channel selection, stability labels or scientific noise definition.
For a new dataset, inspect one representative unit and verify original sample
coordinates and event counts before a large batch.

The CI workflow runs the portable suite on Linux with MATLAB R2024a and Signal
Processing Toolbox using the official [MATLAB Actions](https://github.com/matlab-actions/run-tests).
It does not run GPU/native checks. The repository's Actions page records the
actual outcome of each run.
