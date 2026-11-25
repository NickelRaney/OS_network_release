function mexPath = build_v1_native_mex(verbose)
%BUILD_V1_NATIVE_MEX Expose the standalone C++ engine to MATLAB.
%   Uses the same C++ simulator, L6 and RNG sources as v1_4_standalone.exe.
%   No MATLAB random-number or numerical callbacks occur during simulation.
if nargin < 1, verbose = false; end
validateattributes(verbose, {'logical', 'numeric'}, {'scalar'});
folder = fileparts(mfilename('fullpath'));
compiler = mex.getCompilerConfigurations('C++', 'Selected');
assert(~isempty(compiler), 'V1Native:NoCompiler', ...
    'Select an installed compiler with mex -setup C++.');
args = {'-R2018a', '-O', '-DMINIZ_NO_ARCHIVE_APIS', '-DMINIZ_NO_DEFLATE_APIS'};
if strcmpi(compiler.Manufacturer, 'Microsoft')
    args = [args, {'COMPFLAGS=$COMPFLAGS /std:c++17 /fp:strict /arch:AVX2', ...
        'OPTIMFLAGS=/O2 /Oy- /DNDEBUG'}];
elseif contains(lower(compiler.Manufacturer), 'gnu') || ...
        contains(lower(compiler.Name), 'clang') || ...
        contains(lower(compiler.Name), 'mingw') || contains(lower(compiler.Name), 'gcc')
    args = [args, {'CXXFLAGS=$CXXFLAGS -std=c++17 -fno-fast-math -ffp-contract=off -mavx2', ...
        'CXXOPTIMFLAGS=-O2 -DNDEBUG', 'CFLAGS=$CFLAGS -O2'}];
else
    error('V1Native:Compiler', 'Floating-point options are not configured for %s.', compiler.Name);
end
if verbose, args = [args, {'-v'}]; end
% miniz is needed only because the shared mat5 translation unit also contains
% file I/O. The MEX adapter converts arrays in memory; it performs no MAT I/O.
sources = {fullfile(folder, 'v1_4_native_mex.cpp'), ...
    fullfile(folder, 'v1_simulator.cpp'), fullfile(folder, 'v1_l6.cpp'), ...
    fullfile(folder, 'mat5.cpp'), fullfile(folder, 'third_party', 'miniz', 'miniz_tinfl.c')};
clear v1_4_native_mex
mex(args{:}, ['-I', folder], '-outdir', folder, '-output', 'v1_4_native_mex', sources{:});
mexPath = fullfile(folder, ['v1_4_native_mex.', mexext]);
fprintf('Built %s with %s (native C++ RNG, AVX2, strict floating point).\n', mexPath, compiler.Name);
end
