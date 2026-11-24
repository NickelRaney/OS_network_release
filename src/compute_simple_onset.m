function out = compute_simple_onset(spike_trials, onset_s, cell_indices, varargin)
%COMPUTE_SIMPLE_ONSET Simple-cell onset latency via periodic sliding-window detection.
%   out = compute_simple_onset(spike_trials, onset_s)
%   out = compute_simple_onset(spike_trials, onset_s, cell_indices, Name, Value, ...)
%
% Inputs:
%   spike_trials - nTrials x nPhases cell array of spike matrices (row 1 =
%                  spike count, rows 2:(n+1) = spike times in seconds).
%   onset_s      - Stimulus onset time in seconds (e.g., 1.5).
%   cell_indices - Optional vector of column indices to evaluate (default: all).
%
% Name-Value options:
%   'PeriodMs'     - Grating cycle duration in ms (default: 100 for TF = 10 Hz)
%   'WindowMs'     - Sliding window width in ms (default: 35)
%   'StepMs'       - Bin/step size in ms (default: 5)
%   'MinTrials'    - Minimum trials with >=1 spike in window (default: min(3, nTrials))
%   'MaxLatencyMs' - Censoring cap in ms after onset (default: 500)
%
% Output struct fields:
%   latency_ms       - 1 x nCells minimum latency across phases (ms after onset)
%   phase_latency_ms - nPhases x nCells per-phase latencies (ms after onset)
%   censored         - 1 x nCells logical mask where latency reached MaxLatencyMs
%   median_ms        - Median latency across evaluated cells

assert(iscell(spike_trials) && ~isempty(spike_trials), ...
    'OSNetwork:InvalidInput', 'spike_trials must be a non-empty nTrials x nPhases cell array.');
[n_trials, n_phases] = size(spike_trials);
n_total = size(spike_trials{1,1}, 2);
if nargin < 3 || isempty(cell_indices)
    cell_indices = 1:n_total;
end

p = inputParser;
addParameter(p, 'PeriodMs', 100);
addParameter(p, 'WindowMs', 35);
addParameter(p, 'StepMs', 5);
addParameter(p, 'MinTrials', min(3, n_trials));
addParameter(p, 'MaxLatencyMs', 500);
parse(p, varargin{:});
opt = p.Results;

onset_ms = onset_s * 1000;
n_cells = numel(cell_indices);
phase_lat = opt.MaxLatencyMs * ones(n_phases, n_cells);

for c = 1:n_cells
    cid = cell_indices(c);
    for ph = 1:n_phases
        trial_spikes = cell(n_trials, 1);
        pooled_post = [];
        for tr = 1:n_trials
            S = spike_trials{tr, ph};
            cnt = S(1, cid);
            if cnt > 0
                t_ms = S(2:cnt+1, cid) * 1000;
                t_post = t_ms(t_ms > onset_ms);
            else
                t_post = [];
            end
            trial_spikes{tr} = t_post;
            pooled_post = [pooled_post; t_post(:)]; %#ok<AGROW>
        end
        if isempty(pooled_post)
            continue;
        end
        t_off = find_peak_window(pooled_post, opt.PeriodMs, opt.WindowMs, opt.StepMs);
        while t_off <= opt.MaxLatencyMs
            lo = onset_ms + t_off;
            hi = lo + opt.WindowMs;
            hit_trials = 0;
            hit_spikes = [];
            for tr = 1:n_trials
                ts = trial_spikes{tr};
                in_win = ts(ts > lo & ts < hi);
                if ~isempty(in_win)
                    hit_trials = hit_trials + 1;
                    hit_spikes = [hit_spikes; in_win(:)]; %#ok<AGROW>
                end
            end
            if hit_trials >= opt.MinTrials
                phase_lat(ph, c) = mean(hit_spikes) - onset_ms;
                break;
            end
            t_off = t_off + opt.PeriodMs;
        end
    end
end

out.phase_latency_ms = phase_lat;
out.latency_ms = min(phase_lat, [], 1);
out.censored = out.latency_ms >= opt.MaxLatencyMs;
out.median_ms = median(out.latency_ms);
out.cell_indices = cell_indices;
end

function t_off = find_peak_window(post_spikes_ms, period_ms, window_ms, step_ms)
n_bins = round(period_ms / step_ms);
win_bins = round(window_ms / step_ms);
h = histcounts(mod(post_spikes_ms, period_ms), 0:step_ms:period_ms);
a = conv([h, h], ones(1, win_bins), 'valid');
a = a(1:n_bins);
p = find(a == max(a), 1);
t_off = (mod(p + win_bins - 1, n_bins) - win_bins) * step_ms;
end
