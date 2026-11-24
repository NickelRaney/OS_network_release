function records=source_manifest(root)
%SOURCE_MANIFEST Relative names and hashes of all code and fixed data.
files=[dir(fullfile(root,'src','*.m'));dir(fullfile(root,'config','*.m')); ...
    dir(fullfile(root,'*.m'));dir(fullfile(root,'examples','*.m'));dir(fullfile(root,'data','*.mat'))];
records=repmat(struct('file','','bytes',0,'sha256',''),1,numel(files));
for k=1:numel(files)
    full=fullfile(files(k).folder,files(k).name);
    records(k).file=strrep(full(numel(root)+2:end),'\','/');
    records(k).bytes=files(k).bytes;
    records(k).sha256=file_sha256(full);
end
end
