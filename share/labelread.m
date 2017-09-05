function L = labelread(fname)
fid = fopen(fname,'r');
fgets(fid);
nV  = sscanf(fgets(fid),'%d');
L   = fscanf(fid,'%d %f %f %f %f\n');
L   = reshape(L, 5, nV)';
fclose(fid) ;

