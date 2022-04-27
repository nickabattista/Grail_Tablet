%-------------------------------------------------------------------------------------------------------------------%
%
% IB2d is an Immersed Boundary Code (IB) for solving fully coupled  
% 	fluid-structure interaction models. This version of the code is based off of
%	Peskin's Immersed Boundary Method Paper in Acta Numerica, 2002.
%
% Author: Nicholas A. Battista
% Email:  battistn@tcnj.edu
% Date Created: May 27th, 2015
% Date Modified: April 27th, 2022
% Institution: TCNJ
%
%--------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FUNCTION: imports all Eulerian Data at a single step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fX_Lag,fY_Lag,fLagMag,fLagNorm,fLagTan] = import_Lagrangian_Force_Data_Pendulum(path,numSim)

analysis_path = pwd;

% read in Mag. of Force %
strChoice = 'fX_Lag'; 
fX_Lag =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Tangential Force %
strChoice = 'fY_Lag'; 
fY_Lag = read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Force %
strChoice = 'fMag'; 
fLagMag =  1;%read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Tangential Force %
strChoice = 'fTan'; 
fLagNorm = 1;%read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Normal Force %
strChoice = 'fNorm';
fLagTan =  1;%read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

cd(analysis_path);

clear analysis_path;

clear strChoice first;

