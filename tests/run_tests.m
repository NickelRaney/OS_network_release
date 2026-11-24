function report=run_tests()
%RUN_TESTS Fast data, input and historical-analysis checks; no V1 trial.
root=fileparts(fileparts(mfilename('fullpath')));
report.environment=version;
b=load(fullfile(root,'data','disordered_lgn.mat'));
assert(isequal(size(b.param_v1.cmee),[13689 13689]));
assert(isequal(size(b.param_v1.cmei),[2401 13689]));
assert(isequal(size(b.param_v1.cmie),[13689 2401]));
assert(isequal(size(b.param_v1.cmii),[2401 2401]));
assert(nnz(b.inde(1,:))==841 && nnz(b.indi(1,:))==145);
assert(find(b.kernel3.tker_off~=0,1)-1==16);
assert(find(b.kernel3.tker_on~=0,1)-1==25);
report.historical_onset=reproduce_historical_onset;

% Independently sum AMPA impulses for a hand-specified history. At the
% 12th E and 7th I spikes the original history tests change their gains.
p=b.param_v1;p.duration=0.4;
lgn=zeros(110,25);times=10:10:240;lgn(1,:)=[24 times];
ex=ex_gen2(p,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn);
gainE=[repmat(1.11,1,11),repmat(1.12,1,8),ones(1,5)];
gainI=[repmat(1.15,1,6),repmat(1.05,1,5),ones(1,13)];
time=(1:4000)'*0.1;expected=zeros(4000,2);
for j=1:24
    elapsed=time-times(j);valid=elapsed>=0;
    impulse=(exp(-elapsed(valid)/2)-exp(-elapsed(valid)/0.5))/1.5;
    expected(valid,:)=expected(valid,:)+impulse*[gainE(j) gainI(j)];
end
report.steady_gain_max_error=max(abs(expected-[ex.exe(:,1),ex.exi(:,1)]),[],'all');
assert(report.steady_gain_max_error<1e-12);
% High-cx BG/direct grating starts at 1.15/1.20; later stages are unchanged.
p.lgn_gain_e_initial=1.15;p.lgn_gain_i_initial=1.20;
high=ex_gen2(p,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn);
expected_high=expected;
for j=1:11
    elapsed=time-times(j);valid=elapsed>=0;
    impulse=(exp(-elapsed(valid)/2)-exp(-elapsed(valid)/0.5))/1.5;
    expected_high(valid,1)=expected_high(valid,1)+impulse*(1.15-1.11);
    if j<=6,expected_high(valid,2)=expected_high(valid,2)+impulse*(1.20-1.15);end
end
report.high_cx_steady_gain_max_error=max(abs(expected_high-[high.exe(:,1),high.exi(:,1)]),[],'all');
assert(report.high_cx_steady_gain_max_error<1e-12);
assert(nnz(ex.exe(:,2:end))==0 && nnz(ex.exi(:,2:end))==0);
assert(isequal(full(sum(ex.conn_matrix_e,1))',b.Cat(:,2)));
assert(isequal(full(sum(ex.conn_matrix_i,1))',b.Cat_i(:,1)));

% Onset must leave each LGN at its initial gain until the first spike
% strictly after onset; exercise the complete 20-spike curve and tail.
lgn=zeros(110,42);times=10:10:410;lgn(1,:)=[41 times];
p.duration=0.5;
ex=ex_gen_onset_test(p,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn,0.1);
% The high-cx steady gains above must not silently change the onset curves.
curveE=([1.25 1.28 1.28 1.278 1.275 1.268 1.257 1.240 1.221 1.200 1.177 1.153 1.128 1.104 1.080 1.058 1.038 1.022 1.008 1]-1)*11/25+1;
curveI=([1.28 1.14 1.07 1.0378 1.0112 1 .99 .982 .975 .975 .977 .980 .985 .988 .990 .993 .995 .998 .999 1]-1)*15/28+1;
gains=[repmat([curveE(1);curveI(1)],1,10),[curveE;curveI],ones(2,11)];
time=(1:5000)'*0.1;expected=zeros(5000,2);
for j=1:numel(times)
    elapsed=time-times(j);valid=elapsed>=0;
    impulse=(exp(-elapsed(valid)/2)-exp(-elapsed(valid)/0.5))/1.5;
    expected(valid,:)=expected(valid,:)+impulse*gains(:,j)';
end
report.onset_gain_max_error=max(abs(expected-[ex.exe(:,1),ex.exi(:,1)]),[],'all');
assert(report.onset_gain_max_error<1e-12);
% The Low-LGN model must select its own network, including zero-LGN cells.
cfg=model_config('low_lgn');low=load(fullfile(root,'data',cfg.data_file));
assert(cfg.f1_ee_scale==1.3 && cfg.f1_ei_scale==1.2);
assert(cfg.lgn_gain_e_initial==1.11 && cfg.lgn_gain_i_initial==1.15);
assert(~isequal(low.param_v1.cmee,b.param_v1.cmee));
assert(~isequal(low.LVwire,b.LVwire));
assert(find(low.kernel3.tker_off~=0,1)-1==16 && find(low.kernel3.tker_on~=0,1)-1==25);
lp=low.param_v1;lp.duration=.4;
zero_lgn=zeros(110,25);
low_ex=ex_gen2(lp,low.Cat,low.Cat_i,low.LVwire,low.LVwire_i,zero_lgn);
assert(isequal(full(sum(low_ex.conn_matrix_e,1))',low.Cat(:,2)));
assert(isequal(full(sum(low_ex.conn_matrix_i,1))',low.Cat_i(:,1)));
re=find(low.inde(1,:));ri=size(low.Cat,1)+find(low.indi(1,:));
lc=low.Cat(re,:);counts=arrayfun(@(n)nnz(lc(:,1)==2 & lc(:,2)==n),0:2);
assert(isequal(counts,[74 81 92]));
% A known spike in each complex category must be counted once in that category.
spikes=zeros(2,16090);selected=zeros(1,3);
for n=0:2,selected(n+1)=re(find(lc(:,1)==2 & lc(:,2)==n,1));end
spikes(1,selected)=1;spikes(2,selected)=.2;
summary=spike_summary(spikes,low.Cat,re,ri,[0 .4],'low_lgn');
assert(isequal(summary.cell_count,[841 594 247 74 81 92 145]));
assert(isequal(summary.spike_count,[3 0 3 1 1 1 0]));
report.low_lgn=struct('passed',true,'complex_nLGN',0:2,'complex_counts',counts, ...
    'distinct_network',true,'summary_category_counts',summary.cell_count);
report.passed=true;
fprintf('Input and historical reference checks passed.\n');
end
