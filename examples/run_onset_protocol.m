function out = run_onset_protocol(preset, mode, n_phases, n_trials)
%RUN_ONSET_PROTOCOL Run multi-phase onset trials and compute response latencies.
%   out = run_onset_protocol('disordered_lgn', 'simple', 2, 4)   % Quick preview
%   out = run_onset_protocol('disordered_lgn', 'simple', 10, 4)  % Full simple protocol
%   out = run_onset_protocol('disordered_lgn', 'complex', 20, 3) % Full complex protocol
%
% Each 2.5 s trial (onset at 1.5 s) takes ~1.5--2 min in serial MATLAB.
% Start with small n_phases / n_trials to preview before launching full runs.

if nargin < 1 || isempty(preset), preset = 'disordered_lgn'; end
if nargin < 2 || isempty(mode), mode = 'simple'; end
if nargin < 3 || isempty(n_phases)
    if strcmp(mode, 'simple'), n_phases = 2; else, n_phases = 3; end
end
if nargin < 4 || isempty(n_trials)
    if strcmp(mode, 'simple'), n_trials = 4; else, n_trials = 2; end
end

setup;
base_cfg = model_config(preset);
base_cfg.condition = 'onset';
base_cfg.duration_s = 2.5;
base_cfg.onset_s = 1.5;
base_cfg.analysis_start_s = 1.5;
base_cfg.direction = [1 0];

root = fileparts(fileparts(mfilename('fullpath')));
net = load(fullfile(root, 'data', base_cfg.data_file), 'Cat', 'inde');
roi_e = find(net.inde(1, :));

if strcmp(mode, 'simple')
    phases_pi = (1:n_phases) * (2 / n_phases); % evenly spaced across [0, 2*pi]
    target_cells = roi_e(net.Cat(roi_e, 1) ~= 2);
elseif strcmp(mode, 'complex')
    rng(20261124, 'twister');
    phases_pi = (rand(1, n_phases) - 0.5) * 0.5 + 0.1; % uniform in [-pi/4, pi/4] + 0.1*pi
    target_cells = roi_e(net.Cat(roi_e, 1) == 2);
else
    error('OSNetwork:Mode', 'Mode must be ''simple'' or ''complex''.');
end

spike_trials = cell(n_trials, n_phases);
seed_base = 2026112500;
for ph = 1:n_phases
    for tr = 1:n_trials
        cfg = base_cfg;
        cfg.phase_pi = phases_pi(ph);
        idx = (ph - 1) * n_trials + tr;
        cfg.lgn_seed = seed_base + 2 * idx - 1;
        cfg.v1_seed = seed_base + 2 * idx;
        folder = run_simulation(cfg);
        S = load(fullfile(folder, 'result.mat'), 'spike');
        spike_trials{tr, ph} = S.spike;
    end
end

if strcmp(mode, 'simple')
    out = compute_simple_onset(spike_trials, base_cfg.onset_s, target_cells);
else
    out = compute_complex_onset(spike_trials, base_cfg.onset_s, target_cells);
end
out.preset = base_cfg.preset;
out.mode = mode;
out.phases_pi = phases_pi;
fprintf('%s (%s): median onset latency = %.1f ms (%d cells)\n', ...
    base_cfg.preset, mode, out.median_ms, numel(target_cells));
end
