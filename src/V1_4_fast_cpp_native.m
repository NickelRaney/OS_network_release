function [res, stats] = V1_4_fast_cpp_native(param, init, ex, seedOrOptions)
%V1_4_FAST_CPP_NATIVE Run the standalone C++17/AVX2 engine from MATLAB.
%   res = V1_4_fast_cpp_native(param, init, ex) uses the existing MATLAB input
%   structs and returns VE, VI, HE, HI, L6, amb, lgn, and spike as double
%   arrays. init may be []. No intermediate files or external processes are used.
%
%   res = V1_4_fast_cpp_native(param, init, ex, seed) selects the native C++
%   random stream (xoshiro256** + SplitMix64 + NumPy Ziggurat). seed must be
%   uint64 or an exact nonnegative double <= flintmax (default: 20261005).
%
%   [res, stats] = V1_4_fast_cpp_native(..., options) accepts a scalar struct
%   with fields 'seed' and 'profile' (logical). MATLAB's RNG is neither read
%   nor modified. Because the native C++ RNG differs from MATLAB's 'twister',
%   individual spike trajectories will not match V1_4_fast.m for the same seed.
%
%   Build once before use:
%     setup
%     build_v1_native_mex

if nargin < 4 || isempty(seedOrOptions)
    options = struct();
elseif isstruct(seedOrOptions)
    options = seedOrOptions;
else
    options = struct('seed', seedOrOptions);
end
cppDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'cpp');
if ~strcmp(fileparts(which('v1_4_native_mex')), cppDir)
    addpath(cppDir, '-begin');
end
if exist('v1_4_native_mex', 'file') ~= 3
    error('OSNetwork:CppNotBuilt', ...
        'C++ MEX binary not found. Build it first with: setup; build_v1_native_mex;');
end
[res, stats] = v1_4_native_mex(param, init, ex, options);
end
