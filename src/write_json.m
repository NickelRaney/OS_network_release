function write_json(filename,value)
%WRITE_JSON Save a UTF-8, human-readable run record.
fid=fopen(filename,'w','n','UTF-8');assert(fid>=0,'Cannot write %s',filename);
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(value,PrettyPrint=true));
end
