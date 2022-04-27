%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: function to analyze FORCE (hier_data_IB2d) data from IBAMR
%           simulations / print them to be opened in ParaView or VisIt in
%           VTK format
%
% Author: Nicholas A. Battista
% Date: 8/24/16
% Github: http://github.org/nickabattista
% Institution: UNC-CH
% Lab: Laura Miller Lab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_IBAMR_hier_data_to_VTK()

% Path to Simulation of Interest
%path_name='/Volumes/Marlin2/Trabeculae/Chamber_Tube/APS_Height_Sweep/Re1_Sims/h_pt00/Re1'; 
path_name=pwd;


% Temporal Information
starting_time=20000; % First time step to be included in data analysis
time_step=20000;     % Interval between outputted time steps (this should match the interval in the input2d file)
final_time=600000;   % Final time step to be included in data analysis
timestep=10^(-6);    % Timestep, dt, from the simulation
time=(starting_time:time_step:final_time)*timestep; % Vector of time from simulation


% For Re-dimensionalization (characteristic values)
density=1000; % fluid dynamic density (from IBAMR)
velocity=0.5; % characteristic velocity (m/s^2)
length=0.1;   % characteristic length (m)


% Make Folder To Store .vtk Files
folder_Name = 'force_VTK_Files';
mkdir(folder_Name);

% Print program information to screen
fprintf('\n\n              <--*** PRINTS hier_IB2d_data to VTK format ***--> \n\n');

% Loop over Simulation Time
for k=starting_time:time_step:final_time

    fprintf('->Printing: %d of %d\n',k,final_time);
    
    % Points to FileID for Specific Timestep
    fileID = fopen([path_name '/hier_data_IB2d/F.' num2str(k)]);  %This line opens the file of interest and gives it an FileID #

    % Read in Force Data
    [F_Lag] = read_in_IBAMR_hier_data(fileID);
    
    % Read in Lagrangian Positions
    fileID = fopen([path_name '/hier_data_IB2d/X.' num2str(k)]);
    [lagPts,Npts] = read_in_IBAMR_hier_data(fileID);

    % Compute Normal/Tangential Vectors
    [nX,nY,sqrtN] = give_Me_Lagrangian_Normal_Vectors(Npts,lagPts);
    [tX,tY,~] = give_Me_Lagrangian_Tangent_Vectors(Npts,nX,nY,sqrtN);
    
    % Project Force Data onto Normal / Tangent Vectors
    [F_Tan,F_Normal] = give_Tangent_and_Normal_Force_Projections(Npts,F_Lag,nX,nY,tX,tY);
    
    % Re-dimensionalize the Forces (since the force is the non-dimensional force, or force coefficient)
    F_Tan = -2*F_Tan/(density*velocity^2*length);
    F_Normal = -2*F_Normal/(density*velocity^2*length);
    
    % Compute Colormap Force Magnitude Scalings
    cutoff_Yes = 0; % Don't cutoff top and bottom
    [F_Tan_Mag,F_Normal_Mag] = give_Force_Magnitude_Scalings(Npts,F_Tan,F_Normal,cutoff_Yes);

    % Print Force Data to VTK Format
    print_Forces_to_VTK(folder_Name,num2str(k),lagPts,F_Lag,F_Tan_Mag,F_Normal_Mag)
    
    %Print in MATLAB
    %print_Lagrangian_Mesh_and_Forces(lagPts,F_Tan_Mag,F_Normal_Mag);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    ENDS MAIN PART OF CODE !!!                         %
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: print force data to VTK format
%           -> each lag. pt (x,y) has associated scalar force with it
%           -> force magnitude, normal force, tangential force
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Forces_to_VTK(folder_Name,strNUM,lagPts,F_Lag,F_Tan_Mag,F_Normal_Mag)

 
cd(folder_Name); %change directory to force_VTK_Files folder

    fMagName = ['fMag.' strNUM '.vtk'];
    fNormalName = ['fNorm.' strNUM '.vtk'];
    fTangentName = ['fTan.' strNUM '.vtk'];

    fLagMag = sqrt( F_Lag(:,1).^2 + F_Lag(:,2).^2 ); % Compute magnitude of forces on boundary

    savevtk_points_with_scalar_data( lagPts, fLagMag, fMagName, 'fMag');
    savevtk_points_with_scalar_data( lagPts, F_Normal_Mag, fNormalName, 'fNorm');
    savevtk_points_with_scalar_data( lagPts, F_Tan_Mag, fTangentName, 'fTan');

cd .. % Get out of force_VTK_Files


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints Lagrangian pt data w/ associated scalar to vtk formated file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function savevtk_points_with_scalar_data( X, scalarArray, filename, vectorName)

%X is matrix of size Nx3

N = length( X(:,1) );


%TRY PRINTING THEM AS UNSTRUCTURED_GRID
file = fopen (filename, 'w');
fprintf(file, '# vtk DataFile Version 2.0\n');
fprintf(file, [vectorName '\n']);
fprintf(file, 'ASCII\n');
fprintf(file, 'DATASET UNSTRUCTURED_GRID\n\n');
%
fprintf(file, 'POINTS %i float\n', N);
for i=1:N
    fprintf(file, '%.15e %.15e %.15e\n', X(i,1),X(i,2),0);
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
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: read in data from hier_data from IBAMR
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F_Lag,Npts] = read_in_IBAMR_hier_data(fileID)

%Scans the file for numbers (skipping 'Processor' lines) and puts them in a cell called C
C=textscan(fileID,'%f','HeaderLines',3,'delimiter',' ','commentstyle','P'); 

Forces=C{1}; %Converts cell C to matrix A

F_Lag(:,1) = Forces(1:2:end); % Converts data to x-directed Forces
F_Lag(:,2) = Forces(2:2:end); % Converts data to y-directed Forces
Npts = length(Forces(1:2:end));    % # of Lagrangian Pts.

fclose(fileID); %close the FileID #


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian Derivatives for Normal/Tangential Vector Computation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xL_s,yL_s] = give_Me_Lagrangian_Derivatives(ds,Npts,X,Y)

xL = X;
yL = Y;

xL_s = zeros(Npts,1);
yL_s = zeros(Npts,1);

for i=1:Npts
    if i==1
       xL_s(1) = ( xL(2) - xL(end) ) / (2*ds); 
       yL_s(1) = ( yL(2) - yL(end) ) / (2*ds);
    elseif i<Npts
       xL_s(i) = ( xL(i+1) - xL(i-1) ) / (2*ds); 
       yL_s(i) = ( yL(i+1) - yL(i-1) ) / (2*ds);
    else
       xL_s(i) = ( xL(1) - xL(end-1) ) / (2*ds); 
       yL_s(i) = ( yL(1) - yL(end-1) ) / (2*ds);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian UNIT Normal Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nX,nY,sqrtN] = give_Me_Lagrangian_Normal_Vectors(Npts,lagPts)

% Lagrangian Pts
X = lagPts(:,1); % x-Lagrangian Values
Y = lagPts(:,2); % y-Lagrangian Values

% Compute Lagrangian Spacing
ds = sqrt( ( X(3)-X(4) )^2 + ( Y(3)-Y(4) )^2 );

% Gives Lagrangian Derivatives
[xL_s,yL_s] = give_Me_Lagrangian_Derivatives(ds,Npts,X,Y);

sqrtN = sqrt( (xL_s).^2 + (yL_s).^2 );

nX = ( yL_s ) ./ sqrtN;
nY = ( -xL_s) ./ sqrtN;

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes Lagrangian UNIT Tangent Vectors
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tX,tY,sqrtN] = give_Me_Lagrangian_Tangent_Vectors(Npts,nX,nY,sqrtN)

% Allocate storage
tX = zeros( size(nX) );
tY = zeros( size(nY) );

% Rotate normal vectors to get tangent vectors
ang = -pi/2; % Rotate CW by 90 degrees
for i=1:Npts
    tX(i) = nX(i)*cos(ang) - nY(i)*sin(ang);
    tY(i) = nX(i)*sin(ang) + nY(i)*cos(ang);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: computes force vector projections onto the tangent and normal
%           vectors! 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
function [F_Tan,F_Normal] = give_Tangent_and_Normal_Force_Projections(Npts,F_Lag,nX,nY,tX,tY)

% Allocate Storage
F_Tan = zeros(Npts,2); F_Normal = F_Tan;

Fx = F_Lag(:,1); % Forces in x-direction
Fy = F_Lag(:,2); % Forces in y-direction

for i=1:Npts

    % Compute dot product between force vector and tangent vector
    tanVec_dotProd = ( Fx(i)*tX(i) + Fy(i)*tY(i) ) / sqrt( tX(i)*tX(i) + tY(i)*tY(i) );
    F_Tan(i,1) = tanVec_dotProd * ( tX(i) );
    F_Tan(i,2) = tanVec_dotProd * ( tY(i) );
    
    % Compute dot product between force vector and normal vector
    normalVec_dotProd = ( Fx(i)*nX(i) + Fy(i)*nY(i) ) / sqrt( nX(i)*nX(i) + nY(i)*nY(i) );
    F_Normal(i,1) = normalVec_dotProd * ( nX(i) );
    F_Normal(i,2) = normalVec_dotProd * ( nY(i) );
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: scales the force matrices by desired percentiles of each in magnitude 
%           for colormap scalings 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MagTan,MagNormal] = give_Force_Magnitude_Scalings(Npts,F_Tan,F_Normal,cutoff_Yes)

% Allocate Storage
MagTan = zeros( length( F_Tan(:,1) ),1 );
MagNormal = MagTan;

% Find the magnitude of the force in each 
for i=1:Npts
    MagTan(i,1) = sqrt( F_Tan(i,1)^2  + F_Tan(i,2)^2 );
    MagNormal(i,1)= sqrt( F_Normal(i,1)^2 + F_Normal(i,2)^2);
end

if cutoff_Yes == 1
    
    % Finds Percentiles for Forces for Cutoff Pts.
    prc90_T = prctile(MagTan,90);
    prc90_N = prctile(MagNormal,90);
    prc10_T = prctile(MagTan,10);
    prc10_N = prctile(MagNormal,10);

    % "Cutoff Threshold for the forces" via if-elseif statements by desired percentiles. 
    for i=1:Npts

        mT = MagTan(i);
        mN = MagNormal(i);

        if mT >= prc90_T
            MagTan(i) = prc10_T;
        elseif mT <= prc10_T
            MagTan(i) = prc10_T;
        end

        if mN >= prc90_N
            MagNormal(i) = prc10_N;
        elseif mN <= prc10_N
            MagNormal(i) = prc10_N;
        end

    end
    
end % Ends scaling


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: prints the Lagrangian Mesh with Scalar Value in MATLAB 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function print_Lagrangian_Mesh_and_Forces(lagPts,F_Tan_Mag,F_Normal_Mag)

    subplot(1,2,1)
    interval = [250:7:1500 2600:7:3400];
    scatter(lagPts(interval,1),lagPts(interval,2),55,F_Tan_Mag(interval),'filled');
    title('F_{Tangent}');
    subplot(1,2,2)
    scatter(lagPts(interval,1),lagPts(interval,2),55,F_Normal_Mag(interval),'filled');
    title('F_{Normal}')
    pause(0.5);
    clf;

