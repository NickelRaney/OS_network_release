function conn = build_l4_connectivity(x, y, Cat, model_type)
%BUILD_L4_CONNECTIVITY Construct periodic-boundary L4 connectivity matrices.
%   conn = build_l4_connectivity(x, y, Cat, model_type)
%
% Generates binary cortical connectivity matrices with truncated Gaussian
% distance dependence on a 0.5 x 0.5 deg torus:
%   Presynaptic E: SD = 200/sqrt(2) um (0.10/sqrt(2) deg), cutoff = 250 um (0.125 deg)
%   Presynaptic I: SD = 125/sqrt(2) um (0.0625/sqrt(2) deg), cutoff = 150 um (0.075 deg)
%
% Outputs (stored as presynaptic row x postsynaptic column, matching param_v1):
%   conn.cmee   - Ne x Ne (E -> E)
%   conn.cmie   - Ne x Ni (E -> I)
%   conn.cmei   - Ni x Ne (I -> E)
%   conn.cmii   - Ni x Ni (I -> I)
%   conn.map_ei - 1 x Ni index of nearest E neuron for each I neuron

if nargin < 4 || isempty(model_type), model_type = 'disordered_lgn'; end
if strcmp(model_type, 'high_cx'), model_type = 'disordered_lgn'; end
assert(ismember(model_type, {'disordered_lgn','low_lgn'}), ...
    'OSNetwork:ModelType', 'model_type must be ''disordered_lgn'' or ''low_lgn''.');

sd_e = 0.20 / sqrt(2) / 2;
cutoff_e = 0.25 / 2;
sd_i = 0.125 / sqrt(2) / 2;
cutoff_i = 0.15 / 2;

ne = size(x, 1);
if strcmp(model_type, 'low_lgn')
    peak_ee = 0.104 * ones(ne, 1);
    peak_ee(Cat(:,1) == 2) = 0.156;
else
    peak_ee = 0.097 * ones(ne, 1);
    peak_ee(Cat(:,1) == 5) = 0.096;
    peak_ee(Cat(:,1) == 6) = 0.0835;
    peak_ee(Cat(:,1) == 2 & Cat(:,2) > 9) = 0.1085;
    peak_ee(Cat(:,1) == 2 & Cat(:,2) >= 7 & Cat(:,2) <= 9) = 0.120;
    peak_ee(Cat(:,1) == 2 & Cat(:,2) >= 2 & Cat(:,2) <= 5) = 0.126;
end

% Transpose post x pre -> pre x post to match V1_4_fast convention
conn.cmee = sample_l4_block(x, x, peak_ee, sd_e, cutoff_e)';
conn.cmie = sample_l4_block(x, y, 0.378, sd_e, cutoff_e)';
conn.cmei = sample_l4_block(y, x, 0.80, sd_i, cutoff_i)';
conn.cmii = sample_l4_block(y, y, 0.80, sd_i, cutoff_i)';
conn.map_ei = nearest_e_map(x, y);
end

function wire = sample_l4_block(pre, post, peak, sd, cutoff)
n_pre = size(pre, 1);
n_post = size(post, 1);
if isscalar(peak)
    peak = peak * ones(n_post, 1);
end
wire = false(n_post, n_pre);
pre_x = pre(:, 1)';
pre_y = pre(:, 2)';
cutoff_sq = cutoff^2;
inv_two_sd_sq = 0.5 / (sd^2);

for i = 1:n_post
    dx = mod(post(i, 1) - pre_x, 0.5);
    dx = min(dx, 0.5 - dx);
    dy = mod(post(i, 2) - pre_y, 0.5);
    dy = min(dy, 0.5 - dy);
    dist_sq = dx.^2 + dy.^2;

    prob = peak(i) * exp(-dist_sq * inv_two_sd_sq);
    prob(dist_sq > cutoff_sq | dist_sq < 1e-12) = 0;

    expected = sum(prob);
    trial = rand(1, n_pre) < prob;
    if expected > 0
        for attempt = 1:50
            if abs(sum(trial) - expected) / expected <= 0.10
                break;
            end
            trial = rand(1, n_pre) < prob;
        end
    end
    wire(i, :) = trial;
end
end

function map_ei = nearest_e_map(pre_e, post_i)
n_i = size(post_i, 1);
map_ei = zeros(1, n_i);
pre_x = pre_e(:, 1)';
pre_y = pre_e(:, 2)';
for i = 1:n_i
    dx = mod(post_i(i, 1) - pre_x, 0.5);
    dx = min(dx, 0.5 - dx);
    dy = mod(post_i(i, 2) - pre_y, 0.5);
    dy = min(dy, 0.5 - dy);
    [~, map_ei(i)] = min(dx.^2 + dy.^2);
end
end
