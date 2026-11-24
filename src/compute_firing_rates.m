function rates = compute_firing_rates(spike, window_s)
%COMPUTE_FIRING_RATES Mean firing rate (Hz) per neuron in (window_s(1), window_s(2)].
%   rates = compute_firing_rates(spike, window_s)
%
% Inputs:
%   spike    - Spike matrix where row 1 holds counts and rows 2:(n+1) hold
%              spike times in seconds.
%   window_s - Two-element vector [start_s, end_s] with end_s > start_s.
%
% Output struct fields:
%   all      - 1 x N vector of firing rates (Hz)
%   e        - 1 x 13689 excitatory rates (present when N == 16090)
%   i        - 1 x 2401 inhibitory rates (present when N == 16090)
%   window_s - Analysis window used

validateattributes(window_s, {'numeric'}, {'vector','numel',2,'finite','increasing'});
dur = diff(window_s);
n_cells = size(spike, 2);
fr = zeros(1, n_cells);
for k = 1:n_cells
    n = spike(1, k);
    if n > 0
        t = spike(2:n+1, k);
        fr(k) = sum(t > window_s(1) & t <= window_s(2)) / dur;
    end
end
rates.all = fr;
rates.window_s = window_s;
if n_cells == 16090
    rates.e = fr(1:13689);
    rates.i = fr(13690:end);
end
end
