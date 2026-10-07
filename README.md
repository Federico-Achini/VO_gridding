# VO Gridding

MATLAB implementation of **vector optimization (VO) for NUFFT gridding**.

This repository contains the MATLAB code and MEX implementation used for vector optimization of the **Non-Uniform Fast Fourier Transform (NUFFT)**. The computationally intensive gridding operation is implemented as a compiled MEX function, with support for optimized and parallel execution.

## Features

- MATLAB implementation of vector-optimized NUFFT gridding.
- Native MEX implementation for computationally intensive gridding operations.
- Optional parallel execution through OpenMP.
- Compiler optimization support.
- MATLAB startup script for configuring the project paths.
- Test code included in the repository.
- MIT licensed.

## Repository Structure

```text
VO_gridding/
├── src/                 # MATLAB/C++ source code
├── test/                # Tests and examples
├── grid_par.mexw64      # Precompiled Windows MEX binary
├── mex_build.m          # MEX compilation script
├── startup.m            # MATLAB project startup script
└── LICENSE              # MIT license
```

## Requirements

The project requires:

- MATLAB
- A supported C/C++ compiler for building the MEX function
- OpenMP support if parallel execution is enabled

The repository currently includes a precompiled Windows MEX binary:

```text
grid_par.mexw64
```

The MEX binary is platform-specific. On other operating systems, or when recompilation is required, the MEX function should be built locally.

## Installation

Clone the repository:

```bash
git clone https://github.com/Federico-Achini/VO_gridding.git
cd VO_gridding
```

Open MATLAB and add the repository to your MATLAB path.

The provided `startup.m` script adds both the `src/` and `test/` directories to the MATLAB path:

```matlab
addpath(genpath('src/'))
addpath(genpath('test/'))
```

You can either run these commands manually or execute:

```matlab
startup
```

> **Note:** The current `startup.m` also contains a machine-specific path used during development. If you use the startup script on another machine, remove or update that path as appropriate.

## Building the MEX Function

The repository provides `mex_build.m` for compiling the native gridding implementation.

From MATLAB, run:

```matlab
mex_build
```

The build script compiles:

```text
grid_par.cpp
```

into:

```text
grid_par
```

The current build configuration enables optimization and parallelization:

```matlab
use_optimized = 1;
use_parallel = 1;
```

When parallel compilation is enabled, the build uses OpenMP together with compiler optimizations such as:

```text
-O3
-ffast-math
-march=native
-fopenmp
```

The exact compiler flags may need to be adapted depending on your operating system, compiler, and MATLAB release.

## Usage

The main computational component is the compiled `grid_par` MEX function.

After compiling the MEX function and adding the project directories to the MATLAB path, the gridding routines can be used from MATLAB together with the source code contained in `src/`.

Example project setup:

```matlab
addpath(genpath('src/'))
addpath(genpath('test/'))
```

The `test/` directory contains code demonstrating and validating the implementation.

For details about the expected inputs and outputs of individual functions, refer to the documentation/comments in the corresponding source files.

## Performance

The gridding implementation includes a native MEX backend to reduce the computational overhead associated with performing the gridding operation entirely in MATLAB.

The provided build configuration enables:

- Compiler optimization
- Native architecture optimizations
- Fast floating-point math
- OpenMP-based parallelization

Performance will depend on the MATLAB version, compiler, CPU architecture, number of threads, and problem size.

For reproducible performance comparisons, make sure that the same compiler and optimization settings are used across experiments.

## Testing

Tests are provided in:

```text
test/
```

Make sure the project directories are on the MATLAB path before running the tests:

```matlab
addpath(genpath('src/'))
addpath(genpath('test/'))
```

Then execute the relevant test scripts or functions contained in `test/`.

## Platform Support

The MATLAB source code is intended to be portable, but the compiled MEX binary is platform-dependent.

The repository currently contains:

```text
grid_par.mexw64
```

which is a Windows 64-bit MEX binary.

For Linux or macOS, compile the MEX implementation locally using:

```matlab
mex_build
```

You may need to modify the compiler flags in `mex_build.m` for your platform.

## License

This project is released under the **MIT License**.

See [LICENSE](LICENSE) for the complete license text.

## Author

**Federico Achini**

## Repository

[https://github.com/Federico-Achini/VO_gridding](https://github.com/Federico-Achini/VO_gridding)
