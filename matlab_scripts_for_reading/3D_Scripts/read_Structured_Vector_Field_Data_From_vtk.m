%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: reads in structured point vector field data to vtk formated file 
%           for 3D applications
%
%
% Author: Nicholas A. Battista
% Date: 2/28/17
% Github: http://github.com/nickabattista
% Institution: UNC-CH
% Lab: Laura Miller Lab
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Returns the desired Eulerian Vector Field Data from 3D .vtk
%           file format. 
%
%           NOTE:  (1) User only needs to specify the filename (in string 
%                      format)!
%                  (2) Put file to be analyzed in same folder as this
%                      script
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,V,W,X,Y,Z,Nx,Ny,Nz] = read_Structured_Vector_Field_Data_From_vtk()

print_info();                 % PRINTS SCRIPT INFO

strName = 'visit_ex_db';      % <---File you wish to analyze (w/o '.vtk' tag)

     % % % The function that reads in the Data % % %
[U,V,W,X,Y,Z,Nx,Ny,Nz] = give_Me_The_Data_Please(strName);


     % % % % % DEFINITIONS OF WHAT GETS PASSED BACK % % % % %
     %
     %                U: x-Velocity Component  
     %                V: y-Velocity Component
     %                W: z-Velocity Component
     %           X,Y,Z: X,Y,Z grid values, respectively
     % Nx,Ny,Nz: grid resolution in each direction, respectively



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % ***** END OF WHERE USER SHOULD CHANGE THINGS ***** %
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Reads in the desired Eulerian data from .vtk format 
%           (the actual i/o script)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
function [U,V,W,X,Y,Z,Nx,Ny,Nz] = give_Me_The_Data_Please(strName)        
        
filename = [strName '.vtk'];  % desired file

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

 % Store grid dimensions info
N = sscanf(str,'%*s %f %f %f %*s',3); 
Nx = N(1); Ny = N(2); Nz = N(3);

% bypass lines in header %
str = fgets(fileID);

%X GRID VALUES: get formatting for reading in data from .vtk in fscanf %
strVec = '%f';
for i=2:Nx
    strVec = [strVec ' %f'];
end

% GRID VALUES: read in X vertices %
[X,count] = fscanf(fileID,strVec,Nx);
if count ~= Nx
   error('\nProblem reading in Eulerian Data.'); 
end

% bypass lines in header %
str = fgets(fileID);
str = fgets(fileID);

%Y GRID VALUES: get formatting for reading in data from .vtk in fscanf %
strVec = '%f';
for i=2:Ny
    strVec = [strVec ' %f'];
end

% read in Y vertices %
[Y,count] = fscanf(fileID,strVec,Ny);
if count ~= Ny
   error('\nProblem reading in Eulerian Data.'); 
end

% bypass lines in header %
str = fgets(fileID);
str = fgets(fileID);

%Z GRID VALUES: get formatting for reading in data from .vtk in fscanf %
strVec = '%f';
for i=2:Ny
    strVec = [strVec ' %f'];
end

% read in Z vertices %
[Z,count] = fscanf(fileID,strVec,Nz);
if count ~= Nz
   error('\nProblem reading in Eulerian Data.'); 
end

% AT THIS POINT WE HAVE READ IN THE X,Y,Z grid points! %

% bypass lines in header %
str = fgets(fileID);
str = fgets(fileID);
str = fgets(fileID);


%EULERIAN DATA: get formatting for reading in data from .vtk in fscanf %
strVec = '%f';
for i=2:3*Nx*Ny*Nz
    strVec = [strVec ' %f'];
end

% read in EULERIAN DATA %
[e_Data,count] = fscanf(fileID,strVec,3*Nx*Ny*Nz);
if count ~= 3*Nx*Ny*Nz
   error('\nProblem reading in Eulerian Data.'); 
end


% reshape the matrix into desired data type %
e_Data = reshape(e_Data, 3, count/3); % Reshape (3*Nx*Nx,1) vector to (Nx*Nx,3) matrix
e_Data = e_Data';                     % Store vertices in new matrix

U = e_Data(:,1);       % Store U data
V = e_Data(:,2);       % Store V data
W = e_Data(:,3);       % Store W data

U = reshape(U,Ny,Nx,Nz);  % Reshape (Nx*Ny*Nz,1) matrix to  (Ny,Nx,Nz)
V = reshape(V,Ny,Nx,Nz);  % Reshape (Nx*Ny*Nz,1) matrix to  (Ny,Nx,Nz)
W = reshape(W,Ny,Nx,Nz);  % Reshape (Nx*Ny*Nz,1) matrix to  (Ny,Nx,Nz) 

fclose(fileID);         % Closes the data file.

clear filename fileID str strVec count analysis_path;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: Print Information About What The Script Does
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_info()

fprintf('_________________________________________________________________________\n\n');
fprintf('This script reads in 3D vector field data in .vtk format and...\n');
fprintf('     ...returns (outputs) the following data: \n');
fprintf('          1. U,V,W (the x,y,z components of the velocity)\n');
fprintf('          2. X,Y,Z (the cartesian grid values\n');
fprintf('          3. Nx,Ny,Nz (the grid resolution in each respective direction\n\n');
fprintf('__________________________________________________________________________\n\n');
