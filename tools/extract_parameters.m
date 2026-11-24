function extract_parameters(source_file, destination_file)
%EXTRACT_PARAMETERS Export only the fixed network needed for a release run.
% Run once by a maintainer with access to the original experiment archive.
% No saved external drive or simulation output is reused.
names = {'param','param_v1','init_v1','kernel3','grid', ...
    'Cat','Cat_i','LVwire','LVwire_i','inde','indi'};
s = load(source_file, names{:});
assert(all(isfield(s,names)), 'Snapshot is missing required variables.');
assert(isequal(size(s.param_v1.cmee),[13689 13689]));
assert(isequal(size(s.param_v1.cmei),[2401 13689]));
% Retain exact values/classes, including the original connectivity storage.
save(destination_file, '-struct', 's', '-v7');
check = load(destination_file);
for k = 1:numel(names)
    assert(isequaln(s.(names{k}),check.(names{k})), ...
        'Round-trip mismatch in %s.', names{k});
end
fprintf('Exported and exactly verified %d variables: %s\n',numel(names),destination_file);
for k = 1:numel(names)
    v=s.(names{k}); info=whos('v');
    fprintf('%s: %s %s, %d bytes in memory\n',names{k},class(v),mat2str(size(v)),info.bytes);
end
disp(s.param);
for key={'see','sie','sei','sii','selgn','silgn','dt','duration','onset'}
    fprintf('param_v1.%s = %g\n',key{1},s.param_v1.(key{1}));
end
end
