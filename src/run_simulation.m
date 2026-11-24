function folder=run_simulation(cfg)
%RUN_SIMULATION Readable, serial MATLAB entry point; all seeds are explicit.
root=fileparts(fileparts(mfilename('fullpath')));
validateattributes(cfg.duration_s,{'numeric'},{'scalar','finite','>=',0.4,'<=',10});
assert(abs(cfg.duration_s*10-round(cfg.duration_s*10))<1e-10, ...
    'OSNetwork:Duration','Legacy allocations require multiples of 0.1 s.');
assert(ismember(cfg.condition,{'bg','grating','onset'}));
assert(cfg.onset_s>=0 && cfg.onset_s<cfg.duration_s);
assert(cfg.analysis_start_s>=0 && cfg.analysis_start_s<cfg.duration_s);
assert(isequal(size(cfg.direction),[1 2]) && all(isfinite(cfg.direction)) && norm(cfg.direction)>0);
validateattributes(cfg.kernel_extra_delay_ms,{'numeric'},{'scalar','integer','nonnegative'});
validateattributes(cfg.lgn_seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});
validateattributes(cfg.v1_seed,{'numeric'},{'scalar','integer','>=',0,'<=',2^32-1});

%% Fixed network and explicit run settings
b=load(fullfile(root,'data',cfg.data_file));
param=b.param;
param_v1=b.param_v1;
param.duration=cfg.duration_s;
param_v1.duration=cfg.duration_s;
param_v1.onset=cfg.onset_s;
assert(param.dt==0.1 && param_v1.dt==0.1 && param.speed==0.004);
keys={'f1_ee_scale','f1_ei_scale','lgn_gain_e_initial','lgn_gain_i_initial', ...
    'onset_gain_e_initial','onset_gain_i_initial'};
for k=1:numel(keys)
    validateattributes(cfg.(keys{k}),{'numeric'},{'scalar','positive','finite'});
    param_v1.(keys{k})=cfg.(keys{k});
end

%% Kernel and spatial profile (0.001 degree samples)
kernel3=b.kernel3;
kernel3.tker_off=[zeros(1,cfg.kernel_extra_delay_ms),kernel3.tker_off];
kernel3.tker_on=[zeros(1,cfg.kernel_extra_delay_ms),kernel3.tker_on];
dir=cfg.direction/norm(cfg.direction);
param.LGN_pos=pos_gen(b.grid,dir);
% LGN3's original tail construction has only been audited at this speed.
assert(round(max(param.LGN_pos(:,1))/param.speed)<param.duration*1000);
input=-sin(pi*((1:100000)/200+cfg.phase_pi));
if strcmp(cfg.condition,'bg'), input=zeros(size(input)); end

%% LGN spikes and input to V1
rng(cfg.lgn_seed,'twister');
random_state.lgn_before=rng;
if strcmp(cfg.condition,'onset')
    [out,lgn]=lgnout(param,kernel3,input,cfg.onset_s);
    ex=ex_gen_onset_test(param_v1,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn,cfg.onset_s);
else
    assert(cfg.onset_s==0,'Only onset trials accept a nonzero onset_s.');
    [out,lgn]=lgnout(param,kernel3,input);
    ex=ex_gen2(param_v1,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn);
end
random_state.lgn_after=rng;

%% V1: pure MATLAB, including L6 and ambient stochastic input
rng(cfg.v1_seed,'twister');
random_state.v1_before=rng;
timer=tic;
res=V1_4_fast(param_v1,b.init_v1,ex);
elapsed_s=toc(timer);
random_state.v1_after=rng;
validate_result(res,param_v1);

%% Save to a fresh directory; retain the complete network spike trains
folder=tempname(fullfile(root,'outputs'));
mkdir(folder);
spike=res.spike(1:max(res.spike(1,:))+1,:);
roi_e=find(b.inde(1,:));
roi_i=size(b.Cat,1)+find(b.indi(1,:));
% Legacy trace rows are written at steps 10,20,...; the last allocated row is unused.
rows=1:round(cfg.duration_s*1000);
trace_time_s=rows/1000;
% L6 rows count proposed spikes in blocks starting at step 2,12,... .
l6_block_start_s=(rows-1)/1000+0.0002;
traces.VE=res.VE(rows,roi_e);
traces.VI=res.VI(rows,roi_i-size(b.Cat,1));
traces.HE=res.HE(rows,[roi_e roi_i]);
traces.HI=res.HI(rows,[roi_e roi_i]);
traces.L6=res.L6(rows,[roi_e roi_i]);
traces.ambient=res.amb(rows,roi_e);
traces.LGN=res.lgn(rows,roi_e);
clear res;
summary=spike_summary(spike,b.Cat,roi_e,roi_i,[cfg.analysis_start_s cfg.duration_s],cfg.preset);
meta.matlab_version=version;
meta.computer=computer;
meta.elapsed_v1_seconds=elapsed_s;
meta.source_manifest=source_manifest(root);
meta.note='Fresh MATLAB run; no historical trial rejection and no latency offset.';
save(fullfile(folder,'result.mat'),'cfg','param','param_v1','kernel3','lgn','out', ...
    'random_state','spike','roi_e','roi_i','trace_time_s','l6_block_start_s','traces','summary','meta','-v7.3');
copyfile(fullfile(root,'src'),fullfile(folder,'src'));
copyfile(fullfile(root,'config'),fullfile(folder,'config'));
write_json(fullfile(folder,'run.json'),struct('config',cfg,'summary',summary,'meta',meta));
plot_run(spike,roi_e,roi_i,cfg,summary,fullfile(folder,'raster.png'));
fprintf('Saved %s | %s | %.1f s | V1 %.2f s\n%s\n', ...
    cfg.preset,cfg.condition,cfg.duration_s,elapsed_s,folder);
fprintf('ROI mean E %.6f Hz; I %.6f Hz over (%.3f, %.3f] s\n', ...
    summary.mean_hz(1),summary.mean_hz(end),summary.window_s);
end
