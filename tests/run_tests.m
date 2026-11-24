function report = run_tests()
%RUN_TESTS Fast checks of network data, LGN gains, analysis tools, and scripts.
root = fileparts(fileparts(mfilename('fullpath')));
report.environment = version;

%% 1. Fixed network snapshots & preset aliases
cfg_d = model_config('disordered_lgn');
cfg_h = model_config('high_cx');
assert(isequal(cfg_d, cfg_h), 'high_cx must alias disordered_lgn.');
assert(strcmp(cfg_d.preset, 'disordered_lgn'));

b = load(fullfile(root,'data','disordered_lgn.mat'));
assert(isequal(size(b.param_v1.cmee), [13689 13689]));
assert(isequal(size(b.param_v1.cmei), [2401 13689]));
assert(isequal(size(b.param_v1.cmie), [13689 2401]));
assert(isequal(size(b.param_v1.cmii), [2401 2401]));
assert(nnz(b.inde(1,:)) == 841 && nnz(b.indi(1,:)) == 145);
assert(find(b.kernel3.tker_off ~= 0, 1) - 1 == 16);
assert(find(b.kernel3.tker_on ~= 0, 1) - 1 == 25);
report.historical_onset = reproduce_historical_onset;

%% 2. Steady and onset LGN gain schedules
p = b.param_v1; p.duration = 0.4;
lgn = zeros(110, 25); times = 10:10:240; lgn(1,:) = [24 times];
ex = ex_gen2(p, b.Cat, b.Cat_i, b.LVwire, b.LVwire_i, lgn);
gainE = [repmat(1.11,1,11), repmat(1.12,1,8), ones(1,5)];
gainI = [repmat(1.15,1,6), repmat(1.05,1,5), ones(1,13)];
time = (1:4000)' * 0.1; expected = zeros(4000, 2);
for j = 1:24
    elapsed = time - times(j); valid = elapsed >= 0;
    impulse = (exp(-elapsed(valid)/2) - exp(-elapsed(valid)/0.5)) / 1.5;
    expected(valid,:) = expected(valid,:) + impulse * [gainE(j) gainI(j)];
end
report.steady_gain_max_error = max(abs(expected - [ex.exe(:,1), ex.exi(:,1)]), [], 'all');
assert(report.steady_gain_max_error < 1e-12);

p.lgn_gain_e_initial = 1.15; p.lgn_gain_i_initial = 1.20;
high = ex_gen2(p, b.Cat, b.Cat_i, b.LVwire, b.LVwire_i, lgn);
expected_high = expected;
for j = 1:11
    elapsed = time - times(j); valid = elapsed >= 0;
    impulse = (exp(-elapsed(valid)/2) - exp(-elapsed(valid)/0.5)) / 1.5;
    expected_high(valid,1) = expected_high(valid,1) + impulse * (1.15 - 1.11);
    if j <= 6, expected_high(valid,2) = expected_high(valid,2) + impulse * (1.20 - 1.15); end
end
report.disordered_lgn_steady_gain_max_error = max(abs(expected_high - [high.exe(:,1), high.exi(:,1)]), [], 'all');
assert(report.disordered_lgn_steady_gain_max_error < 1e-12);
assert(nnz(ex.exe(:,2:end)) == 0 && nnz(ex.exi(:,2:end)) == 0);
assert(isequal(full(sum(ex.conn_matrix_e,1))', b.Cat(:,2)));
assert(isequal(full(sum(ex.conn_matrix_i,1))', b.Cat_i(:,1)));

lgn = zeros(110, 42); times = 10:10:410; lgn(1,:) = [41 times];
p.duration = 0.5;
ex = ex_gen_onset_test(p, b.Cat, b.Cat_i, b.LVwire, b.LVwire_i, lgn, 0.1);
curveE = ([1.25 1.28 1.28 1.278 1.275 1.268 1.257 1.240 1.221 1.200 1.177 1.153 1.128 1.104 1.080 1.058 1.038 1.022 1.008 1]-1)*11/25+1;
curveI = ([1.28 1.14 1.07 1.0378 1.0112 1 .99 .982 .975 .975 .977 .980 .985 .988 .990 .993 .995 .998 .999 1]-1)*15/28+1;
gains = [repmat([curveE(1); curveI(1)], 1, 10), [curveE; curveI], ones(2, 11)];
time = (1:5000)' * 0.1; expected = zeros(5000, 2);
for j = 1:numel(times)
    elapsed = time - times(j); valid = elapsed >= 0;
    impulse = (exp(-elapsed(valid)/2) - exp(-elapsed(valid)/0.5)) / 1.5;
    expected(valid,:) = expected(valid,:) + impulse * gains(:,j)';
end
report.onset_gain_max_error = max(abs(expected - [ex.exe(:,1), ex.exi(:,1)]), [], 'all');
assert(report.onset_gain_max_error < 1e-12);

%% 3. Low-LGN snapshot & zero-LGN accounting
cfg = model_config('low_lgn');
low = load(fullfile(root,'data',cfg.data_file));
assert(cfg.f1_ee_scale == 1.3 && cfg.f1_ei_scale == 1.2);
assert(cfg.lgn_gain_e_initial == 1.11 && cfg.lgn_gain_i_initial == 1.15);
assert(~isequal(low.param_v1.cmee, b.param_v1.cmee));
assert(~isequal(low.LVwire, b.LVwire));
assert(find(low.kernel3.tker_off ~= 0, 1) - 1 == 16 && find(low.kernel3.tker_on ~= 0, 1) - 1 == 25);
lp = low.param_v1; lp.duration = 0.4;
low_ex = ex_gen2(lp, low.Cat, low.Cat_i, low.LVwire, low.LVwire_i, zeros(110, 25));
assert(isequal(full(sum(low_ex.conn_matrix_e,1))', low.Cat(:,2)));
assert(isequal(full(sum(low_ex.conn_matrix_i,1))', low.Cat_i(:,1)));
re = find(low.inde(1,:)); ri = size(low.Cat,1) + find(low.indi(1,:));
lc = low.Cat(re,:); counts = arrayfun(@(n) nnz(lc(:,1)==2 & lc(:,2)==n), 0:2);
assert(isequal(counts, [74 81 92]));
spikes = zeros(2, 16090); selected = zeros(1, 3);
for n = 0:2, selected(n+1) = re(find(lc(:,1)==2 & lc(:,2)==n, 1)); end
spikes(1, selected) = 1; spikes(2, selected) = 0.2;
summary = spike_summary(spikes, low.Cat, re, ri, [0 0.4], 'low_lgn');
assert(isequal(summary.cell_count, [841 594 247 74 81 92 145]));
assert(isequal(summary.spike_count, [3 0 3 1 1 1 0]));
report.low_lgn = struct('passed',true,'complex_nLGN',0:2,'complex_counts',counts, ...
    'distinct_network',true,'summary_category_counts',summary.cell_count);

%% 4. Steady-state analysis functions (firing rates, MR, DS, OP)
% Cell 1: purely periodic 10 Hz spikes (MR -> 2); Cell 2: uniform 100 Hz (MR -> 0)
t_periodic = (1.025 : 0.1 : 2.95)';  % 20 spikes over 2 s -> 10 Hz, 1 per cycle
t_uniform  = (1.005 : 0.01 : 2.995)'; % 200 spikes over 2 s -> 100 Hz, flat
S_fwd = zeros(201, 2);
S_fwd(1,1) = numel(t_periodic); S_fwd(2:numel(t_periodic)+1, 1) = t_periodic;
S_fwd(1,2) = numel(t_uniform);  S_fwd(2:numel(t_uniform)+1,  2) = t_uniform;
S_rev = zeros(201, 2);
S_rev(1,1) = 4; S_rev(2:5, 1) = [1.2; 1.7; 2.2; 2.7]; % 2 Hz in opposite direction
S_rev(1,2) = numel(t_uniform); S_rev(2:numel(t_uniform)+1, 2) = t_uniform;

rates_test = compute_firing_rates(S_fwd, [1.0, 3.0]);
assert(abs(rates_test.all(1) - 10) < 1e-12 && abs(rates_test.all(2) - 100) < 1e-12);

mr_test = compute_modulation_ratio(S_fwd, [1.0, 3.0], 10, S_rev);
assert(abs(mr_test.mr_fwd(1) - 2.0) < 1e-12);
assert(abs(mr_test.mr_fwd(2)) < 1e-12);
assert(abs(mr_test.mr_pref(1) - 2.0) < 1e-12 && abs(mr_test.f0_pref(1) - 10) < 1e-12);

fr_f = 20 * ones(1, 13689); fr_r = 10 * ones(1, 13689); fr_o = 5 * ones(1, 13689);
sel_test = compute_selectivity(fr_f, fr_r, fr_o, b.Cat, find(b.inde(1,:)), 'disordered_lgn');
assert(all(abs(sel_test.median_ds - 2.0) < 1e-12));
assert(all(abs(sel_test.median_op - 0.25) < 1e-12));

%% 5. Onset latency estimators (simple & complex)
% Synthetic simple cell firing at 1.540 s, 1.640 s, ... across all 4 trials (40 ms latency)
simp_trials = cell(4, 2);
for tr = 1:4
    S1 = zeros(11, 1); ts1 = (1.540 : 0.1 : 2.440)';
    S1(1,1) = numel(ts1); S1(2:numel(ts1)+1, 1) = ts1;
    S2 = zeros(11, 1); ts2 = (1.570 : 0.1 : 2.470)';
    S2(1,1) = numel(ts2); S2(2:numel(ts2)+1, 1) = ts2;
    simp_trials{tr, 1} = S1;
    simp_trials{tr, 2} = S2;
end
simp_out = compute_simple_onset(simp_trials, 1.5, 1);
assert(abs(simp_out.latency_ms(1) - 40) < 1e-9);

% Synthetic complex cell: 0 Hz before 1.550 s, 100 Hz after 1.550 s -> T50 ~ 50 ms
cx_trials = cell(1, 10);
ts_cx = (1.5505 : 0.01 : 2.495)';
for tr = 1:10
    Sc = zeros(numel(ts_cx)+1, 1);
    Sc(1,1) = numel(ts_cx); Sc(2:end, 1) = ts_cx;
    cx_trials{tr} = Sc;
end
cx_out = compute_complex_onset(cx_trials, 1.5, 1);
assert(cx_out.latency_ms(1) >= 40 && cx_out.latency_ms(1) <= 60);

%% 6. Network construction script smoke checks
[lgn_grid, gx, gy] = build_lgn_grid();
assert(isequal(size(lgn_grid), [2 110]));
assert(isequal(size(gx), [13689 2]) && isequal(size(gy), [2401 2]));
tpl = build_lgn_templates();
assert(isequal(size(tpl.temp), [3 4 2]));
rng(123, 'twister');
[lv_sub, cat_sub, lvi_sub, cati_sub] = assign_lgn_connections(gx(1:50,:), gy(1:20,:), tpl, 'disordered_lgn');
assert(isequal(size(lv_sub), [50 14]) && isequal(size(cat_sub), [50 3]));
assert(isequal(size(lvi_sub), [20 7]) && isequal(size(cati_sub), [20 1]));
assert(max(lv_sub(:)) <= 110 && max(lvi_sub(:)) <= 110);
conn_sub = build_l4_connectivity(gx(1:50,:), gy(1:20,:), cat_sub, 'disordered_lgn');
assert(isequal(size(conn_sub.cmee), [50 50]) && isequal(size(conn_sub.cmie), [50 20]));
assert(isequal(size(conn_sub.cmei), [20 50]) && isequal(size(conn_sub.cmii), [20 20]));
assert(isequal(size(conn_sub.map_ei), [1 20]));

report.analysis_and_scripts_passed = true;
report.passed = true;
fprintf('All network, input, analysis, and script checks passed.\n');
end
