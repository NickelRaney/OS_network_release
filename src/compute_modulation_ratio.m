function res = compute_modulation_ratio(spike_fwd, window_s, tf_hz, spike_rev)
%COMPUTE_MODULATION_RATIO Compute modulation ratio MR = F1/F0 from drifting-grating spikes.
%   res = compute_modulation_ratio(spike_fwd, window_s)
%   res = compute_modulation_ratio(spike_fwd, window_s, tf_hz)
%   res = compute_modulation_ratio(spike_fwd, window_s, tf_hz, spike_rev)
%
% Uses 10 ms histogram bins folded over integer 1-second epochs and evaluates
% MR = 2 * |F(tf_hz)| / |F(0)| via FFT (default tf_hz = 10 Hz).
% When both opposite optimal directions (spike_fwd and spike_rev) are
% provided, res.mr_pref selects the modulation ratio from each neuron's
% preferred direction (higher mean rate F0).

if nargin < 3 || isempty(tf_hz), tf_hz = 10; end
validateattributes(window_s, {'numeric'}, {'vector','numel',2,'finite','increasing'});
validateattributes(tf_hz, {'numeric'}, {'scalar','integer','>=',1,'<=',49});

dur_s = diff(window_s);
n_sec = round(dur_s);
assert(abs(dur_s - n_sec) < 1e-9 && n_sec >= 1, ...
    'OSNetwork:ModulationWindow', ...
    'Analysis window duration must be a positive integer number of seconds.');

[mr_fwd, f0_fwd] = single_direction_mr(spike_fwd, window_s, n_sec, tf_hz);
res.mr = mr_fwd;
res.f0 = f0_fwd;
res.window_s = window_s;
res.tf_hz = tf_hz;

if nargin >= 4 && ~isempty(spike_rev)
    assert(size(spike_rev, 2) == size(spike_fwd, 2), ...
        'OSNetwork:SizeMismatch', 'spike_fwd and spike_rev must have the same number of columns.');
    [mr_rev, f0_rev] = single_direction_mr(spike_rev, window_s, n_sec, tf_hz);
    use_fwd = f0_fwd >= f0_rev;
    mr_pref = mr_rev;
    mr_pref(use_fwd) = mr_fwd(use_fwd);
    res.mr_fwd = mr_fwd;
    res.mr_rev = mr_rev;
    res.f0_fwd = f0_fwd;
    res.f0_rev = f0_rev;
    res.mr_pref = mr_pref;
    res.f0_pref = max(f0_fwd, f0_rev);
end
end

function [mr, f0] = single_direction_mr(spike, window_s, n_sec, tf_hz)
edges_ms = (window_s(1)*1000) : 10 : (window_s(2)*1000);
n_cells = size(spike, 2);
mr = nan(1, n_cells);
f0 = zeros(1, n_cells);
harm_idx = tf_hz + 1;
for k = 1:n_cells
    n = spike(1, k);
    if n == 0, continue; end
    t_ms = spike(2:n+1, k) * 1000;
    counts = histcounts(t_ms, edges_ms);
    folded = sum(reshape(counts, 100, n_sec), 2);
    total = sum(folded);
    f0(k) = total / n_sec;
    if total > 0
        fc = fft(folded);
        mr(k) = 2 * abs(fc(harm_idx)) / abs(fc(1));
    end
end
end
