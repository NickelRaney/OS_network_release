function hex=file_sha256(filename)
%FILE_SHA256 Stream binary bytes using MATLAB's bundled Java runtime.
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(filename,'rb');assert(fid>=0,'Cannot open %s',filename);
cleanup=onCleanup(@() fclose(fid));
while ~feof(fid)
    bytes=fread(fid,1024*1024,'*uint8');
    md.update(typecast(bytes,'int8'));
end
hex=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
