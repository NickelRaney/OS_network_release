function cfg = model_config(preset)
%MODEL_CONFIG Return default simulation configuration for a model preset.
%
% Supported presets:
%   'disordered_lgn' (default; alias: 'high_cx') - Disordered-LGN model
%   'low_lgn'                                    - Low-LGN model
%
% Stimulus direction convention for the central ROI (vertical-preferring):
%   [1 0] or [-1 0]  - vertical drifting grating (optimal orientation)
%   [0 1] or [0 -1]  - horizontal drifting grating (orthogonal orientation)

if nargin==0 || isempty(preset)
    preset = 'disordered_lgn';
end
preset = char(preset);
if strcmp(preset,'high_cx')
    preset = 'disordered_lgn';
end

cfg.preset = preset;
cfg.data_file = 'disordered_lgn.mat';
switch cfg.preset
    case 'disordered_lgn'
        cfg.f1_ee_scale = 2.25;
        cfg.f1_ei_scale = 1.5;
        cfg.lgn_gain_e_initial = 1.15;
        cfg.lgn_gain_i_initial = 1.20;
    case 'low_lgn'
        cfg.data_file = 'low_lgn.mat';
        cfg.f1_ee_scale = 1.3;
        cfg.f1_ei_scale = 1.2;
        cfg.lgn_gain_e_initial = 1.11;
        cfg.lgn_gain_i_initial = 1.15;
    otherwise
        error('OSNetwork:Preset', ...
            'Unknown preset "%s". Use ''disordered_lgn'' (alias ''high_cx'') or ''low_lgn''.', ...
            cfg.preset);
end

% Condition-specific onset initial gains (used by ex_gen_onset_test)
cfg.onset_gain_e_initial = 1.11;
cfg.onset_gain_i_initial = 1.15;

% Keep saved OFF/ON kernel padding (16/25 ms) with no extra shift by default
cfg.kernel_extra_delay_ms = 0;

cfg.condition = 'grating'; % 'bg', 'grating', or 'onset'
cfg.duration_s = 0.4;
cfg.onset_s = 0;
cfg.direction = [1 0];     % [1 0] / [-1 0]: optimal (vertical); [0 1]: orthogonal (horizontal)
cfg.phase_pi = 0.1;

% Fixed spatial frequency (2.5 cycles/deg) and speed (4 deg/s -> TF = 10 Hz)
cfg.backend = 'matlab';    % 'matlab' (default) or 'cpp_native' (requires cpp/build_v1_native_mex)
cfg.lgn_seed = 2026112301;
cfg.v1_seed = 2026112302;
cfg.analysis_start_s = 0;
end
