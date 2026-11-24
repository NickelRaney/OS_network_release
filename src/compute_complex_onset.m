function out = compute_complex_onset(spike_phase_groups, onset_s, cell_indices, varargin)
%COMPUTE_COMPLEX_ONSET Complex-cell onset latency via Gaussian-smoothed firing rate.
%   out = compute_complex_onset(spike_phase_groups, onset_s)
%   out = compute_complex_onset(spike_phase_groups, onset_s, cell_indices, Name, Value, ...)
%
% Inputs:
%   spike_phase_groups - Either:
%                        (a) 1 x nGroups cell array, where each entry is a
%                            cell array of spike matrices for that phase group, or
%                        (b) a single cell array of spike matrices (1 phase group).
%   onset_s            - Stimulus onset time in seconds (e.g., 1.5).
%   cell_indices       - Optional vector of column indices to evaluate (default: all).
%
% Name-Value options:
%   'ThresholdFraction' - Fraction a in f_bg + a*(f_ss - f_bg) (default: 0.50)
%   'KernelSDMs'        - Gaussian smoothing SD in ms (default: 5)
%   'SustainMs'         - Subsequent duration in ms that must stay above threshold (default: 25)
%   'MaxSearchMs'       - Maximum post-onset latency searched in ms (default: 300)
%   'PreWindowMs'       - Pre-onset baseline window in ms (default: 1000)
%   'PostWindowMs'      - Post-onset window in ms (default: 1000)
%   'SteadyStartMs'     - Post-onset start of steady-state window in ms (default: 200)
%
% Output struct fields:
%   latency_ms       - 1 x nCells fastest latency across phase groups (MaxSearchMs+1 if none)
%   phase_latency_ms - nGroups x nCells latencies per phase group
%   winning_phase    - 1 x nCells index of winning phase group
%   bg_hz            - 1 x nCells background firing rate (winning phase)
%   steady_hz        - 1 x nCells steady-state firing rate (winning phase)
%   no_crossing      - 1 x nCells logical mask where threshold was never crossed
%   median_ms        - Median latency across evaluated cells

assert(iscell(spike_phase_groups) && ~isempty(spike_phase_groups), ...
    'OSNetwork:InvalidInput', 'spike_phase_groups must be a non-empty cell array.');
if ~iscell(spike_phase_groups{1})
    spike_phase_groups = {spike_phase_groups};
end

first_mat = spike_phase_groups{1}{1};
n_total = size(first_mat, 2);
if nargin < 3 || isempty(cell_indices)
    cell_indices = 1:n_total;
end

p = inputParser;
addParameter(p, 'ThresholdFraction', 0.50);
addParameter(p, 'KernelSDMs', 5);
addParameter(p, 'SustainMs', 25);
addParameter(p, 'MaxSearchMs', 300);
addParameter(p, 'PreWindowMs', 1000);
addParameter(p, 'PostWindowMs', 1000);
addParameter(p, 'SteadyStartMs', 200);
parse(p, varargin{:});
opt = p.Results;

n_groups = numel(spike_phase_groups);
n_cells = numel(cell_indices);
pre_bins = opt.PreWindowMs;
post_bins = opt.PostWindowMs;
edges_s = (onset_s - pre_bins/1000) : 0.001 : (onset_s + post_bins/1000);

half_w = round(2 * opt.KernelSDMs);
if opt.KernelSDMs == 5, half_w = 10; end % 21-point kernel matching paper/legacy
k_grid = -half_w:half_w;
g_kernel = exp(-(k_grid.^2) / (2 * opt.KernelSDMs^2));
g_kernel = g_kernel / sum(g_kernel);

no_cross_val = opt.MaxSearchMs + 1;
phase_lat = no_cross_val * ones(n_groups, n_cells);
bg_all = zeros(n_groups, n_cells);
ss_all = zeros(n_groups, n_cells);

for g = 1:n_groups
    trials = spike_phase_groups{g}(:);
    n_tr = numel(trials);
    rate_scale = 1000 / n_tr;
    for c = 1:n_cells
        cid = cell_indices(c);
        pooled = [];
        for tr = 1:n_tr
            S = trials{tr};
            cnt = S(1, cid);
            if cnt > 0
                pooled = [pooled; S(2:cnt+1, cid)]; %#ok<AGROW>
            end
        end
        counts = histcounts(pooled, edges_s);
        inst_fr = conv(counts * rate_scale, g_kernel, 'same');
        bg = mean(inst_fr(1:pre_bins));
        ss = mean(inst_fr((pre_bins + opt.SteadyStartMs):end));
        bg_all(g, c) = bg;
        ss_all(g, c) = ss;
        thresh = bg + opt.ThresholdFraction * (ss - bg);
        max_t = min(opt.MaxSearchMs, post_bins - opt.SustainMs);
        for t = 1:max_t
            if min(inst_fr((pre_bins + t):(pre_bins + t + opt.SustainMs))) > thresh
                phase_lat(g, c) = t;
                break;
            end
        end
    end
end

[best_lat, win_idx] = min(phase_lat, [], 1);
bg_win = zeros(1, n_cells);
ss_win = zeros(1, n_cells);
for c = 1:n_cells
    bg_win(c) = bg_all(win_idx(c), c);
    ss_win(c) = ss_all(win_idx(c), c);
end

out.latency_ms = best_lat;
out.phase_latency_ms = phase_lat;
out.winning_phase = win_idx;
out.bg_hz = bg_win;
out.steady_hz = ss_win;
out.no_crossing = best_lat >= no_cross_val;
out.median_ms = median(best_lat);
out.cell_indices = cell_indices;
end
