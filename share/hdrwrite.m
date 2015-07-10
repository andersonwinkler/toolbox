function hdrwrite(hdr,filehdr,endianness)
% Overwrite the header of a NIFTI file.
% 
% hdrwrite(hdr,filename,endianness)
% 
% - filename   : File to be read. It can be a
%                .hdr or a .nii file.
% - hdr        : Struct containing the header.
%                See the 'hdrread' command.
% - endianness : Define the endianness:
%                'l': little endian
%                'b': big endian
% 
% _____________________________________
% Anderson M. Winkler
% FMRIB / Univ. of Oxford
% Nov/2012
% http://brainder.org

% Write down (overwrite!)
fid=fopen(filehdr,'r+',endianness);
fwrite(fid,hdr.sizeof_hdr,'int');
fwrite(fid,hdr.data_type,'char');
fwrite(fid,hdr.db_name,'char');
fwrite(fid,hdr.extents,'int');
fwrite(fid,hdr.session_error,'short');
fwrite(fid,hdr.regular,'char');
fwrite(fid,hdr.dim_info,'char');
fwrite(fid,hdr.dim,'short');
fwrite(fid,hdr.intent_p1,'float');
fwrite(fid,hdr.intent_p2,'float');
fwrite(fid,hdr.intent_p3,'float');
fwrite(fid,hdr.intent_code,'short');
fwrite(fid,hdr.datatype,'short');
fwrite(fid,hdr.bitpix,'short');
fwrite(fid,hdr.slice_start,'short');
fwrite(fid,hdr.pixdim,'float');
fwrite(fid,hdr.vox_offset,'float');
fwrite(fid,hdr.scl_slope,'float');
fwrite(fid,hdr.scl_inter,'float');
fwrite(fid,hdr.slice_end,'short');
fwrite(fid,hdr.slice_code,'char');
fwrite(fid,hdr.xyzt_units,'char');
fwrite(fid,hdr.cal_max,'float');
fwrite(fid,hdr.cal_min,'float');
fwrite(fid,hdr.slice_duration,'float');
fwrite(fid,hdr.toffset,'float');
fwrite(fid,hdr.glmax,'int');
fwrite(fid,hdr.glmin,'int');
fwrite(fid,hdr.descrip,'char');
fwrite(fid,hdr.aux_file,'char');
fwrite(fid,hdr.qform_code,'short');
fwrite(fid,hdr.sform_code,'short');
fwrite(fid,hdr.quatern_b,'float');
fwrite(fid,hdr.quatern_c,'float');
fwrite(fid,hdr.quatern_d,'float');
fwrite(fid,hdr.qoffset_x,'float');
fwrite(fid,hdr.qoffset_y,'float');
fwrite(fid,hdr.qoffset_z,'float');
fwrite(fid,hdr.srow_x,'float');
fwrite(fid,hdr.srow_y,'float');
fwrite(fid,hdr.srow_z,'float');
fwrite(fid,hdr.intent_name,'char');
fwrite(fid,hdr.magic,'char');
fclose(fid);