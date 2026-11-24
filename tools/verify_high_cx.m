function report=verify_high_cx(archive_root)
%VERIFY_HIGH_CX Compare with the original expressions plus confirmed constants.
% This verifies the release implementation, not the historical paper seeds.
root=fileparts(fileparts(mfilename('fullpath')));
scratch=tempname(fullfile(root,'outputs'));mkdir(scratch);
text=fileread(fullfile(archive_root,'V1_4_fast.m'));
text=strrep(text,'function [res]=V1_4_fast(','function [res]=corrected_reference(');
assert(contains(text,')*1.3/3+0.92;') && contains(text,')*1.2/1.5+0.95;'));
text=strrep(text,')*1.3/3+0.92;',')*2.25/3+0.92;');
text=strrep(text,')*1.2/1.5+0.95;',')*1.5/1.5+0.95;');
fid=fopen(fullfile(scratch,'corrected_reference.m'),'w','n','UTF-8');
fprintf(fid,'%s',text);fclose(fid);
addpath(scratch);cleanup=onCleanup(@() rmpath(scratch));
b=load(fullfile(root,'data','disordered_lgn.mat'));
p=b.param_v1;p.duration=.4;p.onset=0;
cfg=model_config('high_cx');
keys={'f1_ee_scale','f1_ei_scale','lgn_gain_e_initial','lgn_gain_i_initial'};
for k=1:numel(keys),p.(keys{k})=cfg.(keys{k});end
b.param.duration=.4;b.param.LGN_pos=pos_gen(b.grid,[1 0]);
rng(cfg.lgn_seed,'twister');[~,lgn]=lgnout(b.param,b.kernel3,-sin(pi*((1:100000)/200+cfg.phase_pi)));
ex=ex_gen2(p,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn);
rng(cfg.v1_seed,'twister');ref=corrected_reference(p,b.init_v1,ex);state=rng;
rng(cfg.v1_seed,'twister');candidate=V1_4_fast(p,b.init_v1,ex);
assert(isequal(state,rng));
names=fieldnames(ref);
for k=1:numel(names),assert(isequaln(ref.(names{k}),candidate.(names{k})),'Mismatch %s',names{k});end
report.passed=true;report.exact_output_fields=names;report.exact_final_rng=true;
report.duration_s=.4;report.f1_scales=[2.25 1.5];report.steady_lgn_initial=[1.15 1.20];
report.note='Original MATLAB reference with only the two author-confirmed F1 literals changed.';
report.matlab_version=version;
report.source_manifest=source_manifest(root);
write_json(fullfile(root,'validation','high_cx_comparison.json'),report);
fprintf('High-cx corrected-reference comparison: all eight output fields and RNG exactly equal.\n');
end
