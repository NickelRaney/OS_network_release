function folder = run_demo(preset)
%RUN_DEMO Short simulation of the full 13,689 E / 2,401 I network.
% This is a dependency/replay demonstration, not a steady-state estimate.
if nargin==0, preset='high_cx'; end
setup;
cfg=model_config(preset);
cfg.condition='grating';
cfg.duration_s=0.4;
cfg.onset_s=0;
cfg.analysis_start_s=0;
cfg.lgn_seed=2026112301;
cfg.v1_seed=2026112302;
folder=run_simulation(cfg);
end
