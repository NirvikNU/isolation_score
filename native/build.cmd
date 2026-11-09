@echo off
rem Build from an x64 Visual Studio Developer Command Prompt.
rem Set CUDA_PATH to a toolkit with CUDA 11 cuBLAS and NVRTC (validated: 11.2).
if not defined CUDA_PATH (echo Set CUDA_PATH first. & exit /b 1)
cl /nologo /O2 /openmp /arch:AVX2 /fp:precise /std:c++17 /EHsc /DBOUND_DIMS=32 /DEXACT_OWN_AP /DCONTIGUOUS_LEAVES isolation_knn.cpp /Fe:isolation_knn.exe
if errorlevel 1 exit /b 1
cl /nologo /O2 /fp:precise /std:c++17 /EHsc /I"%CUDA_PATH%\include" isolation_is.cpp /Fe:isolation_is.exe /link /LIBPATH:"%CUDA_PATH%\lib\x64" cuda.lib nvrtc.lib cublas.lib
exit /b %errorlevel%
