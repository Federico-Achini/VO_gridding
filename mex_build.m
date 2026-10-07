% --------------------
% |  compile_grid.m   |
% --------------------
%
% This script compiles the spline gridding MEX function using the
% configured optimization and parallelization options.
%
% The source file and output name are defined at the beginning of the
% script, while the compiler flags are selected according to the
% optimization and parallelization settings.
%
% Author: Federico Achini
% Date:   August 2025
%
% ------------------------------------------------------------------------

%% Compilation options

use_optimized = 1;
use_parallel = 1;


%% Source and output configuration

path = '';
name = 'grid_par';

filename = strcat(path, name, '.cpp');


%% Locate source file

fprintf('=== Compiling spline gridding MEX (Windows/MSVC) ===\n');

srcRoot = pwd;
srcFile = fullfile(srcRoot, filename);

if ~exist(srcFile, 'file')

    error('Source file not found: %s', srcFile);

end


%% Configure compiler arguments

fprintf('Compiling with MSVC...\n');

if use_optimized

    if use_parallel

        mexArgs = { ...

            fullfile(srcRoot, filename), ...

            '-R2018a', ...

            'CXXFLAGS=$CXXFLAGS  -ffast-math -march=native -fopenmp -O3 -save-temps -v -O3', ...

            'LDFLAGS=$LDFLAGS -fopenmp', ...

            '-output', ...

            fullfile(srcRoot, name) ...

        };

    else

        mexArgs = { ...

            fullfile(srcRoot, filename), ...

            '-R2018a', ...

            'CFLAGS=$CFLAGS -O3 -ffast-math -march=native', ...

            'LDFLAGS=$LDFLAGS', ...

            '-output', ...

            fullfile(srcRoot, name) ...

        };

    end

else

    mexArgs = { ...

        fullfile(srcRoot, filename), ...

        '-R2018a', ...

        'COMPFLAGS=$COMPFLAGS /O2 /arch:AVX2 /MT', ...

        'LINKFLAGS=$LINKFLAGS /NODEFAULTLIB:msvcrt', ...

        '-output', ...

        fullfile(srcRoot, name) ...

    };

end


%% Compile MEX function

mex(mexArgs{:});