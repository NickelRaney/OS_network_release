function folders = run_conditions(preset)
%RUN_CONDITIONS Run single BG, optimal grating, and onset trials in sequence.
%   folders = run_conditions()                uses 'disordered_lgn'.
%   folders = run_conditions('disordered_lgn') (alias: 'high_cx').
%   folders = run_conditions('low_lgn')        uses the Low-LGN model.
if nargin==0, preset = 'disordered_lgn'; end
cfg = model_config(preset);
cfg.duration_s = 2.5;
cfg.lgn_seed = 2026112311;
cfg.v1_seed = 2026112312;
folders = cell(1, 3);
conditions = {'bg','grating','onset'};
for k = 1:3
    cfg.condition = conditions{k};
    cfg.onset_s = 0;
    cfg.analysis_start_s = 1;
    if strcmp(cfg.condition,'onset')
        cfg.onset_s = 1.5;
        cfg.analysis_start_s = 1.5;
    end
    folders{k} = run_simulation(cfg);
end
end
