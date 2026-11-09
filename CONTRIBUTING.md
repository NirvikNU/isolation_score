# Contributing

Please open an issue with a small synthetic reproduction, MATLAB/platform
version, options, population counts, status and error identifier. Do not post
recordings or identifying metadata that you cannot share publicly.

Public APIs live directly in `+isolation/`. Internal numerical functions live
in `+isolation/private/`; optional native sources live in `native/`. Keep
acquisition-specific loaders outside the core. The examples illustrate the
input contract without requiring private files.

Before submitting a change:

1. Run `runtests('tests/portable')` and `assertSuccess` on the results.
2. For changes to native code, build it and run `tests/native` on supported
   hardware. Check both mixed and double IS precision.
3. Compare scientific changes with an independent exhaustive reference on
   small inputs, including ties and rare populations. A matching wrapper
   test alone does not establish mathematical correctness.
4. Update documentation and add a focused regression for changed behavior.
5. If changing the metric definition, event selection or preprocessing,
   document the difference and revise `MetricVersion` deliberately. Do not
   silently reinterpret saved scores.

The repository uses the MIT License. Include appropriate attribution for new
third-party code and ensure it can be distributed under compatible terms.
