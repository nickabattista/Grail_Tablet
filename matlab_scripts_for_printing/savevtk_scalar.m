%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints 3D scalar matrix to VTK formated file
%
% Author: Nicholas A. Battista
% Created: 08/24/16
% Modified: 04/27/22
% Github: http://github.org/nickabattista
% Lab: TCNJ Bioinspiration Lab
% Institution: TCNJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_scalar(array, filename, colorMap, dx, dy, time)

% array is matrix of size 'nx x ny x nz' containing scalar data on
%              computational grid
% filename:    What you are saving the VTK file as (string)
% colorMap:    What you are naming the data you're printing (string)
% dx,dy:       Grid spacing (resolution)
% time: current time in simulation

[ny, nx, nz] = size(array);
fid = fopen(filename, 'wt');
fprintf(fid, '# vtk DataFile Version 2.0\n');
fprintf(fid, 'Comment goes here\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, '\n');
fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
%
fprintf(fid, 'FIELD FieldData 1\n');
fprintf(fid, 'TIME 1 1 double\n');
fprintf(fid, '%.8f\n',time);
%
fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', ny, nx, nz);
fprintf(fid, '\n');
fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
%fprintf(fid, 'SPACING   1.000   1.000   1.000\n'); if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
fprintf(fid, ['SPACING   ' num2str(dx) ' '   num2str(dy) '   1.000\n']);
fprintf(fid, '\n');
fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
fprintf(fid, ['SCALARS ' colorMap ' double\n']);
fprintf(fid, 'LOOKUP_TABLE default\n');
fprintf(fid, '\n');
for a=1:nz
    for b=1:nx
        for c=1:ny
            fprintf(fid, '%d ', array(c,b,a));
        end
        fprintf(fid, '\n');
    end
end
fclose(fid);
    
return