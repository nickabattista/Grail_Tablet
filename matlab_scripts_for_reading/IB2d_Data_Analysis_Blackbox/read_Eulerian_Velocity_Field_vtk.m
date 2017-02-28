%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  nick.battista@unc.edu
% Date Created: May 27th, 2015
% Institution: UNC-CH
%
% This code is capable of creating Lagrangian Structures using:
% 	1. Springs
% 	2. Beams (*torsional springs)
% 	3. Target Points
%	4. Muscle-Model (combined Force-Length-Velocity model, "HIll+(Length-Tension)")
%
% One is able to update those Lagrangian Structure parameters, e.g., spring constants, resting lengths, etc
% 
% There are a number of built in Examples, mostly used for teaching purposes. 
% 
% If you would like us to add a specific muscle model, please let Nick (nick.battista@unc.edu) know.
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the velocity data field from .vtk format
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,V] = read_Eulerian_Velocity_Field_vtk(path,simNums)

analysis_path = pwd; % Store path to analysis folder!

cd(path);

filename = ['u.' num2str(simNums) '.vtk'];  % desired lagPts.xxxx.vtk file

fileID = fopen(filename);
if ( fileID== -1 )
    error('\nCANNOT OPEN THE FILE!');
end

str = fgets(fileID); %-1 if eof
if ~strcmp( str(3:5),'vtk');
    error('\nNot in proper VTK format');
end

% read in the header info %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);

% Store grid info
N = sscanf(str,'%*s %f %*f %*s',1); 
Nx = N(1);

% bypass lines in header %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);


% get formatting for reading in data from .vtk in fscanf %
strVec = '%f';
for i=2:3*Nx
    strVec = [strVec ' %f'];
end

% read in the vertices %
[e_Data,count] = fscanf(fileID,strVec,3*Nx*Nx);
if count ~= 3*Nx*Nx
   error('Problem reading in Eulerian Data.'); 
end

% reshape the matrix into desired data type %
e_Data = reshape(e_Data, 3, count/3); % Reshape (3*Nx*Nx,1) vector to (Nx*Nx,3) matrix
e_Data = e_Data';                     % Store vertices in new matrix

U = e_Data(:,1);       % Store U data
V = e_Data(:,2);       % Store V data

U = reshape(U,Nx,Nx)';  % Reshape (Nx*Nx,1) matrix to  (Nx,Nx)
V = reshape(V,Nx,Nx)';  % Reshape (Nx*Nx,1) matrix to  (Nx,Nx)
 
fclose(fileID);         % Closes the data file.

cd ..;                  % Change directory back to ../viz_IB2d/ directory

cd(analysis_path);      % Change directory back to Data Analysis Folder

clear filename fileID str strVec count analysis_path e_Data Nx;

