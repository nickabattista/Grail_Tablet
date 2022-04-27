%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints unstructured point data w/ associated scalar value to vtk 
%           formated file
%
% Author: Nicholas A. Battista
% Created: 08/24/16
% Modified: 04/27/22
% Github: http://github.org/nickabattista
% Lab: TCNJ Bioinspiration Lab
% Institution: TCNJ
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_points_with_scalar_data( X, scalarArray, filename, vectorName, time)

% X is matrix of size Nx3 
%              Col 1: x-data
%              Col 2: y-data
%              Col 3: z-data
% scalarArray: Scalar array you are assigning to each point
% filename:    What you are saving the VTK file as (string)
% vectorName:  What you are naming the data you're printing (string)
% time: current time in simulation

N = length( X(:,1) );


%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'FIELD FieldData 1\n');
fprintf(file, 'TIME 1 1 double\n');
fprintf(file, '%.8f\n',time);
%
fprintf(file, 'POINTS %i float\n', N);
for i=1:N
    fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),X(i,3));
end
fprintf(file,'\n');
%
fprintf(file, 'POINT_DATA   %d\n', N);
fprintf(file, ['SCALARS ' vectorName ' double\n']);
fprintf(file, 'LOOKUP_TABLE default\n');
fprintf(file, '\n');
    for i=1:N
        fprintf(file, '%d ', scalarArray(i,1));
        fprintf(file, '\n');
    end

fclose(file);
    
    