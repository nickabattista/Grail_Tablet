%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints 3D scalar matrix to VTK formated file
%
%
% Author: Nicholas A. Battista
% Date: 8/24/16
% Github: http://github.org/nickabattista
% Institution: UNC-CH
% Lab: Laura Miller Lab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_scalar(array, filename, colorMap,dx,dy)

% array is matrix of size 'nx x ny x nz' containing scalar data on
%              computational grid
% filename:    What you are saving the VTK file as (string)
% colorMap:    What you are naming the data you're printing (string)
% dx,dy:       Grid spacing (resolution)

    [nx, ny, nz] = size(array);
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    fprintf(fid, ['SPACING   ' num2str(dx) ' '   num2str(dy) '   1.000\n']);
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, ['SCALARS ' colorMap ' double\n']);
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    for a=1:nz
        for b=1:ny
            for c=1:nx
                fprintf(fid, '%d ', array(c,b,a));
            end
            fprintf(fid, '\n');
        end
    end
    fclose(fid);
return