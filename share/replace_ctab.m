function replace_ctab(oldannotfile,ctabfile,newannotfile)
% Replaces the colortable of a FreeSurfer annotation file. The
% result is saved in a new annotation file.
% 
% Usage:
% replace_ctab(oldannotfile,ctabfile,newannotfile)
% 
% Inputs:
% oldannotfile : Annotation file which colortable is going to be replaced.
% ctabfile     : Customised color table to be used in the new file.
% newannotfile : New annotation file to be created.
% 
% Before running, be sure that ${FREESURFER_HOME}/matlab is in your
% OCTAVE/MATLAB path.
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Jul/2011

% Read the annotation file which colortable is going to be changed
[vtx,lab,oldctab] = read_annotation(oldannotfile);

% Read the colortable that is going to be used
fid = fopen(ctabfile);
C = textscan(fid,'%d %s %d %d %d %d\n');
fclose(fid);

% Reformat the new colortable
newctab.numEntries = numel(C{1});
newctab.orig_tab = ctabfile;
newctab.struct_names = C{2};
newctab.table = double([C{3} C{4} C{5} C{6}]);

% Add the column for the structure ID
newctab.table(:,5) = newctab.table(:,1:4) * [1 2^8 2^16 2^24]';

% Initialize a new label var
labn = zeros(size(lab));

% For each structure (search by name), replace the ID in the label variable
for s = 1:newctab.numEntries,
    
    % Structure OLD index & ID
    sidx = find(strcmp(newctab.struct_names{s},oldctab.struct_names));
    osid = oldctab.table(sidx,5);
    
    % Structure NEW ID
    nsid = newctab.table(s,5);
    
    % Replace OLD by NEW IDs in the label var
    labn(lab == osid) = nsid;
end

% Save as the new annotation file
write_annotation(newannotfile,vtx,labn,newctab);
