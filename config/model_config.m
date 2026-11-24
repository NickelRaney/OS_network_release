function cfg = model_config(preset)
%MODEL_CONFIG Explicit historical constants; see docs/PROVENANCE.md.
if nargin==0, preset='high_cx'; end
cfg.preset=char(preset);
cfg.data_file='disordered_lgn.mat';
switch cfg.preset
    case 'low_lgn'
        % Separate Low-LGN network and constants from parameters.txt.
        cfg.data_file='low_lgn.mat';
        cfg.f1_ee_scale=1.3;
        cfg.f1_ei_scale=1.2;
        cfg.lgn_gain_e_initial=1.11;
        cfg.lgn_gain_i_initial=1.15;
    case 'historical_disk'
        % High-cx network with disk constants; NOT the Low-LGN network.
        cfg.f1_ee_scale=1.3;
        cfg.f1_ei_scale=1.2;
        cfg.lgn_gain_e_initial=1.11;
        cfg.lgn_gain_i_initial=1.15;
    case 'high_cx'
        % Values listed under "New model" in the original parameters.txt.
        cfg.f1_ee_scale=2.25;
        cfg.f1_ei_scale=1.5;
        cfg.lgn_gain_e_initial=1.15;
        cfg.lgn_gain_i_initial=1.20;
    otherwise
        error('OSNetwork:Preset','Unknown preset: %s',cfg.preset);
end
% Preserve the separately audited historical onset curves, not the steady
% history rule. These are condition-specific settings of this same model.
cfg.onset_gain_e_initial=1.11;
cfg.onset_gain_i_initial=1.15;
% Keep the saved OFF/ON padding (16/25 ms). No implicit latency shift.
cfg.kernel_extra_delay_ms=0;
cfg.condition='grating'; % bg, grating or onset
cfg.duration_s=0.4;
cfg.onset_s=0;
cfg.direction=[1 0];
cfg.phase_pi=0.1;
% Fixed historical SF=2.5 cycles/degree and speed=4 degrees/second (TF=10 Hz).
cfg.lgn_seed=2026112301;
cfg.v1_seed=2026112302;
cfg.analysis_start_s=0;
end
