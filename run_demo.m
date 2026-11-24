function folder = run_demo(preset)
%RUN_DEMO Short 0.4 s simulation of the full 13,689 E / 2,401 I V1 network.
%   folder = run_demo()                runs the Disordered-LGN model.
%   folder = run_demo('disordered_lgn') runs the Disordered-LGN model (alias: 'high_cx').
%   folder = run_demo('low_lgn')        runs the Low-LGN model.
%
% This short trial demonstrates execution and output generation; use longer
% durations (e.g., 3 s with analysis_start_s = 1) for steady-state estimates.
if nargin==0, preset = 'disordered_lgn'; end
setup;
cfg = model_config(preset);
cfg.condition = 'grating';
cfg.duration_s = 0.4;
cfg.onset_s = 0;
cfg.analysis_start_s = 0;
cfg.lgn_seed = 2026112301;
cfg.v1_seed = 2026112302;
folder = run_simulation(cfg);
end
