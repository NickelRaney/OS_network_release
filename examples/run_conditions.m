function folders=run_conditions(preset)
%RUN_CONDITIONS One fresh BG, grating and onset trial, without trial selection.
if nargin==0, preset='high_cx'; end
cfg=model_config(preset);
cfg.duration_s=2.5;
cfg.lgn_seed=2026112311;
cfg.v1_seed=2026112312;
folders=cell(1,3);
conditions={'bg','grating','onset'};
for k=1:3
    cfg.condition=conditions{k};
    cfg.onset_s=0;
    cfg.analysis_start_s=1;
    if strcmp(cfg.condition,'onset')
        cfg.onset_s=1.5;
        cfg.analysis_start_s=1.5;
    end
    folders{k}=run_simulation(cfg);
end
end
