function report=verify_against_archive(archive_root,preset)
%VERIFY_AGAINST_ARCHIVE Maintainer-only comparison to untouched MATLAB files.
% Call after setup. The original experiment directory is read only.
root=fileparts(fileparts(mfilename('fullpath')));
if nargin<2,preset='historical_disk';end
assert(ismember(preset,{'historical_disk','low_lgn'}), ...
    'Use verify_high_cx for the author-corrected high-cx reference.');
cfg=model_config(preset);
scratch=tempname(fullfile(root,'outputs'));mkdir(scratch);
names={'V1_4_fast','ex_gen2','ex_gen_onset_test'};
for k=1:numel(names)
    text=fileread(fullfile(archive_root,[names{k} '.m']));
    % Rename only the public function, preserving every expression/statement.
    text=regexprep(text,'(?m)^function \[res\]=\w+\(', ...
        ['function [res]=original_' names{k} '('],'once');
    fid=fopen(fullfile(scratch,['original_' names{k} '.m']),'w','n','UTF-8');
    fprintf(fid,'%s',text);fclose(fid);
end
addpath(scratch);cleanup=onCleanup(@() rmpath(scratch));
b=load(fullfile(root,'data',cfg.data_file));
b.param.duration=0.4;b.param_v1.duration=0.4;b.param_v1.onset=0;
b.param.LGN_pos=pos_gen(b.grid,[1 0]);
rng(2026112301,'twister');
[~,lgn]=lgnout(b.param,b.kernel3,-sin(pi*((1:100000)/200+0.1)));
ex=ex_gen2(b.param_v1,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn);
ref=original_ex_gen2(b.param_v1,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn);
assert(isequal(ex,ref));report.exact_steady_input=true;clear ref;
onset=ex_gen_onset_test(b.param_v1,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn,.1);
ref=original_ex_gen_onset_test(b.param_v1,b.Cat,b.Cat_i,b.LVwire,b.LVwire_i,lgn,.1);
assert(isequal(onset,ref));report.exact_onset_input=true;clear onset ref;

rng(2026112302,'twister');timer=tic;
ref=original_V1_4_fast(b.param_v1,b.init_v1,ex);reference_seconds=toc(timer);
ref_rng=rng;
fprintf('Original MATLAB completed in %.2f s.\n',reference_seconds);
% Test the explicit disk preset, including both new F1 fields.
b.param_v1.f1_ee_scale=cfg.f1_ee_scale;b.param_v1.f1_ei_scale=cfg.f1_ei_scale;
rng(2026112302,'twister');timer=tic;
candidate=V1_4_fast(b.param_v1,b.init_v1,ex);candidate_seconds=toc(timer);
candidate_rng=rng;
assert(isequal(ref_rng,candidate_rng));
names=fieldnames(ref);
for k=1:numel(names)
    assert(isequaln(ref.(names{k}),candidate.(names{k})),'Mismatch: %s',names{k});
end
validate_result(candidate,b.param_v1);
report.passed=true;report.exact_output_fields=names;
report.preset=preset;report.data_file=cfg.data_file;
report.exact_final_rng=true;report.duration_s=.4;
report.lgn_seed=2026112301;report.v1_seed=2026112302;
report.original_seconds=reference_seconds;report.release_seconds=candidate_seconds;
report.matlab_version=version;report.computer=computer;
report.reference_simulator_sha256=file_sha256(fullfile(archive_root,'V1_4_fast.m'));
report.release_simulator_sha256=file_sha256(fullfile(root,'src','V1_4_fast.m'));
destination='archive_comparison.json';
if strcmp(preset,'low_lgn'),destination='low_lgn_archive_comparison.json';end
write_json(fullfile(root,'validation',destination),report);
fprintf('Original/release: all eight output fields and final RNG are exactly equal.\n');
end
