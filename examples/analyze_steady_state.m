function metrics = analyze_steady_state(preset, duration_s)
%ANALYZE_STEADY_STATE Run BG, optimal (+/-), and orthogonal gratings and compute metrics.
%   metrics = analyze_steady_state()                   % Disordered-LGN, 3 s trials
%   metrics = analyze_steady_state('low_lgn', 3.0)     % Low-LGN, 3 s trials
%
% Computes firing rates (after a 1 s startup transient), modulation ratio
% (MR = F1/F0), directional selectivity (DS), and orthogonal/preferred (OP)
% ratio for the central ROI (vertical-preferring hypercolumn).

if nargin < 1 || isempty(preset), preset = 'disordered_lgn'; end
if nargin < 2 || isempty(duration_s), duration_s = 3.0; end
assert(duration_s >= 2.0 && abs(duration_s - round(duration_s)) < 1e-9, ...
    'OSNetwork:Duration', 'Use an integer duration >= 2 s for MR FFT analysis.');

setup;
base_cfg = model_config(preset);
base_cfg.duration_s = duration_s;
base_cfg.onset_s = 0;
base_cfg.analysis_start_s = 1.0;

% Conditions: BG, optimal forward [1 0], optimal reverse [-1 0], orthogonal [0 1]
specs = { ...
    'bg',      'bg',      [1 0],  2026112401, 2026112402; ...
    'opt_fwd', 'grating', [1 0],  2026112411, 2026112412; ...
    'opt_rev', 'grating', [-1 0], 2026112421, 2026112422; ...
    'ortho',   'grating', [0 1],  2026112431, 2026112432};

folders = struct();
spikes = struct();
rates = struct();
for k = 1:size(specs, 1)
    name = specs{k, 1};
    cfg = base_cfg;
    cfg.condition = specs{k, 2};
    cfg.direction = specs{k, 3};
    cfg.lgn_seed = specs{k, 4};
    cfg.v1_seed = specs{k, 5};
    folders.(name) = run_simulation(cfg);
    S = load(fullfile(folders.(name), 'result.mat'), 'spike', 'roi_e', 'roi_i');
    spikes.(name) = S.spike;
    rates.(name) = compute_firing_rates(S.spike, [1.0, duration_s]);
    roi_e = S.roi_e;
end

root = fileparts(fileparts(mfilename('fullpath')));
net = load(fullfile(root, 'data', base_cfg.data_file), 'Cat');

metrics.preset = base_cfg.preset;
metrics.folders = folders;
metrics.window_s = [1.0, duration_s];
metrics.rates = rates;
metrics.selectivity = compute_selectivity( ...
    rates.opt_fwd.e, rates.opt_rev.e, rates.ortho.e, net.Cat, roi_e, base_cfg.preset);
metrics.modulation_ratio = compute_modulation_ratio( ...
    spikes.opt_fwd, [1.0, duration_s], 10, spikes.opt_rev);

simple_roi = roi_e(net.Cat(roi_e, 1) ~= 2);
complex_roi = roi_e(net.Cat(roi_e, 1) == 2);
fprintf('\n=== Steady-state summary (%s, window [1, %.0f] s) ===\n', ...
    base_cfg.preset, duration_s);
fprintf('Simple median MR = %.3f | Complex median MR = %.3f\n', ...
    median(metrics.modulation_ratio.mr_pref(simple_roi), 'omitnan'), ...
    median(metrics.modulation_ratio.mr_pref(complex_roi), 'omitnan'));
fprintf('Simple median DS = %.3f | Complex median DS = %.3f\n', ...
    metrics.selectivity.median_ds(2), metrics.selectivity.median_ds(3));
fprintf('Simple median OP = %.3f | Complex median OP = %.3f\n', ...
    metrics.selectivity.median_op(2), metrics.selectivity.median_op(3));
end
