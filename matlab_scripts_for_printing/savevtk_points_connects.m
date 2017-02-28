%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints unstructured grid (CONNECTED point) data to VTK formated file
%
%
% Author: Nicholas A. Battista
% Date: 8/24/16
% Github: http://github.org/nickabattista
% Institution: UNC-CH
% Lab: Laura Miller Lab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_points_connects(X, filename, vectorName,connectsMat)

% X is matrix of size Nx3 
%              Col 1: x-data
%              Col 2: y-data
%              Col 3: z-data
% filename:    What you are saving the VTK file as (string)
% vectorName:  What you are naming the data you're printing (string)
% connectsMat: Describes connections between nodes 

N = length( X(:,1) );
Nc = length( connectsMat(:,1) );

% PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'POINTS %i float\n', N);
for i=1:N
    fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),X(i,3));
end
fprintf(file,'\n');
%
fprintf(file,'CELLS %i %i\n',Nc,3*Nc); %First: # of "Cells", Second: Total # of info inputed following
for s=1:Nc
    fprintf(file,'%i %i %i\n',2, connectsMat(s,1), connectsMat(s,2) );
end
fprintf(file,'\n');
%
fprintf(file,'CELL_TYPES %i\n',Nc); % N = # of "Cells"
for i=1:Nc
   fprintf(file,'3 '); 
end
fprintf(file,'\n');
fclose(file);