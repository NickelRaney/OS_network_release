function root = setup()
%SETUP Add only this release's code to the MATLAB path.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'),fullfile(root,'config'), ...
    fullfile(root,'examples'),fullfile(root,'tests'));
assert(~isempty(which('exprnd')), ...
    'OSNetwork:Toolbox', 'NLIF requires Statistics and Machine Learning Toolbox (exprnd).');
assert(strcmp(which('V1_4_fast'),fullfile(root,'src','V1_4_fast.m')), ...
    'OSNetwork:Shadowing','Another V1_4_fast shadows this release. Change to the release root.');
end
