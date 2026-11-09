# Changelog

## 1.0.0 — 2026-11-09

- Public MATLAB package with explicit raw-recording and prepared-waveform APIs.
- Stable interval and acquisition-boundary validation independent of any file
  format, behavioral event schema or session table.
- Full-population IS and paper-normalized FN/FP; all-AP SNR statistics.
- Optional CPU/CUDA acceleration preserving all waveform coordinates.
- Sequential batch jobs with lazy loading, resume and progress reporting.
- Independently selected unit diagnostics and result summaries.
- Synthetic examples, exhaustive reference tests, native checks and CI.

The public adapter preserves the established calculation implementation.
Packaging changes relocate native resources and make input/plot labels generic;
they do not change metric formulas or preprocessing defaults.
