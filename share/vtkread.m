function [Vtx,Fac] = vtkread(filename)
% This is not supposed to read any other VTK file, but only
% those such as created by FSL-FIRST (i.e. POLYDATA).
%
% Usage:
% [vtx,fac] = vtkread(filename)
%
% _____________________________________
% Anderson M. Winkler
% Yale University / Institute of Living
% Nov/2009

fid = fopen(filename);
T = textscan(fid,'%s','CommentStyle','#','Delimiter','*'); % fake delimiter
fclose(fid);

for t = 1:10,
    if strcmp(T{1,1}{t,1}(1:5),'POINT'),
        nVtx = textscan(T{1,1}{t,1},'%s');
        nVtx = str2num(nVtx{1,1}{2,1});
        Vtx = zeros(nVtx,3);
        for v = 1:nVtx,
           tmpvar = textscan(T{1,1}{v+t,1},'%n');
           Vtx(v,:) = tmpvar{1}';
        end
        blah = t+nVtx+1;
        nFac = numel(T{1,1})-blah;
        Fac = zeros(nFac,3);
        for f  = 1:nFac,
            tmpvar = textscan(T{1,1}{f+blah,1},'%n');
            Fac(f,:) = tmpvar{1}(2:4)';
        end
    end
end
Fac = Fac+1;