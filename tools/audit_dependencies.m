function report=audit_dependencies()
%AUDIT_DEPENDENCIES Check MATLAB's static dependency closure of public entries.
root=fileparts(fileparts(mfilename('fullpath')));
entries={fullfile(root,'run_demo.m'),fullfile(root,'examples','run_conditions.m'), ...
    fullfile(root,'examples','reproduce_historical_onset.m'),fullfile(root,'tests','run_tests.m')};
[files,products]=matlab.codetools.requiredFilesAndProducts(entries);
relative=cell(size(files));
for k=1:numel(files)
    if startsWith(files{k},[root filesep])
        relative{k}=strrep(files{k}(numel(root)+2:end),'\','/');
    else
        assert(startsWith(files{k},matlabroot), ...
            'OSNetwork:ExternalDependency','Unexpected external dependency: %s',files{k});
        relative{k}=strrep(files{k}(numel(matlabroot)+2:end),'\','/');
    end
end
report.files=relative;report.products=products;report.passed=true;
report.note='Static source closure; dynamic data paths are covered by executed demos and packaged data.';
write_json(fullfile(root,'validation','dependencies.json'),report);
disp({products.Name});
end
