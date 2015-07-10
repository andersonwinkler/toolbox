function hdr = hdrread(filehdr,endianness)
% Read the header of a NIFTI file.
% 
% hdr = hdrread(filename,endianness)
% 
% - filename   : File to be read. It can be a
%                .hdr or a .nii file.
% - endianness : Define the endianness:
%                'l': little endian
%                'b': big endian
% - hdr        : Struct containing the header.
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Nov/2012
% http://brainder.org

% Read the header
fid=fopen(filehdr,'r',endianness);
hdr.sizeof_hdr     = fread(fid,1,'int');
hdr.data_type      = fread(fid,10,'char');
hdr.db_name        = fread(fid,18,'char');
hdr.extents        = fread(fid,1,'int');
hdr.session_error  = fread(fid,1,'short');
hdr.regular        = fread(fid,1,'char');
hdr.dim_info       = fread(fid,1,'char');
hdr.dim            = fread(fid,8,'short');
hdr.intent_p1      = fread(fid,1,'float');
hdr.intent_p2      = fread(fid,1,'float');
hdr.intent_p3      = fread(fid,1,'float');
hdr.intent_code    = fread(fid,1,'short');
hdr.datatype       = fread(fid,1,'short');
hdr.bitpix         = fread(fid,1,'short');
hdr.slice_start    = fread(fid,1,'short');
hdr.pixdim         = fread(fid,8,'float');
hdr.vox_offset     = fread(fid,1,'float');
hdr.scl_slope      = fread(fid,1,'float');
hdr.scl_inter      = fread(fid,1,'float');
hdr.slice_end      = fread(fid,1,'short');
hdr.slice_code     = fread(fid,1,'char');
hdr.xyzt_units     = fread(fid,1,'char');
hdr.cal_max        = fread(fid,1,'float');
hdr.cal_min        = fread(fid,1,'float');
hdr.slice_duration = fread(fid,1,'float');
hdr.toffset        = fread(fid,1,'float');
hdr.glmax          = fread(fid,1,'int');
hdr.glmin          = fread(fid,1,'int');
hdr.descrip        = fread(fid,80,'char');
hdr.aux_file       = fread(fid,24,'char');
hdr.qform_code     = fread(fid,1,'short');
hdr.sform_code     = fread(fid,1,'short');
hdr.quatern_b      = fread(fid,1,'float');
hdr.quatern_c      = fread(fid,1,'float');
hdr.quatern_d      = fread(fid,1,'float');
hdr.qoffset_x      = fread(fid,1,'float');
hdr.qoffset_y      = fread(fid,1,'float');
hdr.qoffset_z      = fread(fid,1,'float');
hdr.srow_x         = fread(fid,4,'float');
hdr.srow_y         = fread(fid,4,'float');
hdr.srow_z         = fread(fid,4,'float');
hdr.intent_name    = fread(fid,16,'char');
hdr.magic          = fread(fid,4,'char');
fclose(fid);

% Some fields to char
hdr.db_name     = char(hdr.db_name');
hdr.descrip     = char(hdr.descrip');
hdr.aux_file    = char(hdr.aux_file');
hdr.intent_name = char(hdr.intent_name');