%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints structured point vector data to vtk formated file
%
% Author: Nicholas A. Battista
% Created: 08/24/16
% Modified: 04/27/22
% Github: http://github.org/nickabattista
% Lab: TCNJ Bioinspiration Lab
% Institution: TCNJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_vector(X, Y, filename, vectorName, dx, dy, time)

% X,Y is matrix of size Nx3 containing X,Y-direction vector components
%              Col 1: x-data
%              Col 2: y-data
%              Col 3: z-data
% filename:    What you are saving the VTK file as (string)
% vectorName:  What you are naming the data you're printing (string)
% dx,dy:       Grid spacing (resolution)
% time: current time in simulation

if (size(X) ~= size(Y))
    fprint('Error: velocity arrays of unequal size\n'); return;
end
[nx, ny, nz] = size(X); 
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
fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
fprintf(fid, '\n');
fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
%fprintf(fid, 'SPACING   1.000   1.000   1.000\n'); if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly]
fprintf(fid, ['SPACING   ' num2str(dx) ' '  num2str(dy) '   1.000\n']);
fprintf(fid, '\n');
fprintf(fid, 'POINT_DATA   %d\n', nx*ny);
fprintf(fid, ['VECTORS ' vectorName ' double\n']);
fprintf(fid, '\n');
for a=1:nz
    for b=1:ny
        for c=1:nx
            fprintf(fid, '%f ', X(c,b,1));
            fprintf(fid, '%f ', Y(c,b,1));
            fprintf(fid, '%f ', 0);
        end
        fprintf(fid, '\n');
    end

end
fclose(fid);
return

