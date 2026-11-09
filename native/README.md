# Optional native acceleration

The optional native backend accelerates scoring for large recordings using
CPU neighbor searches and GPU affinity calculations. The MATLAB backend can
be used without these executables.

## Supported configuration

- 64-bit Windows and an AVX2-capable CPU.
- Visual Studio 2022 C++ build tools, including a Windows SDK and OpenMP.
- NVIDIA GPU with compute capability 6.0 or newer, an appropriate driver, and
  the CUDA toolkit. The tested toolkit is **CUDA 11.2**; other versions are not
  currently validated.
- Exactly 144 waveform columns, K=31, and more than 31 pooled events.
- Sufficient RAM, temporary disk, and GPU memory for the selected populations.

MATLAB Parallel Computing Toolbox is not required: MATLAB launches standalone
executables. Source is supplied; generated binaries and vendor runtime DLLs
are deliberately excluded from Git.

## Build

Open an **x64 Native Tools Command Prompt for Visual Studio 2022**, then:

```bat
cd /d C:\path\to\isolation_score\native
set "CUDA_PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2"
build.cmd
```

This produces `isolation_knn.exe` and `isolation_is.exe` next to the sources.
Keep `isolation_kernels.cu` alongside the executables: IS compiles these kernels
at runtime through NVRTC. Set the same toolkit root in MATLAB:

```matlab
options = struct('Backend','native', ...
    'CUDARoot','C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2');
[result, diagnostics] = isolation.score(recording, options);
```

Alternatively set `CUDA_PATH` before starting MATLAB. Native tests use that
environment variable to locate the toolkit:

```matlab
setenv('CUDA_PATH','C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2');
results = runtests('tests/native');
assertSuccess(results);
```

Missing tools are reported as an error. Do not download arbitrary replacement
executables. Linux/macOS native builds are not currently supported by the
MATLAB adapter, even if parts of the C++ sources can be ported.

## Implementation

`isolation_knn.cpp` implements CPU FN/FP, with conservative search bounds and
full-coordinate FP64 final distances. `isolation_is.cpp` uses CUDA/cuBLAS and
`isolation_kernels.cu` for all-pair AP affinities. The default GPU path uses
mixed precision; `ISPrecision='double'` requests FP64 affinities.

The build flags in `build.cmd` are part of the tested configuration. In
particular, the 32-coordinate orthonormal bound, exact target search, contiguous
leaves and precise floating-point setting should not be removed casually.

Native scratch data consists of row-major little-endian FP64 waveforms and a
projection matrix, written in a new temporary directory. These are internal
transport files, not a public acquisition format. The adapter validates output
population counts, score normalization and mean-affinity consistency, then
removes its temporary directory.
