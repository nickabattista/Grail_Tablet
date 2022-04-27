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

function [fX_Lag,fY_Lag] = import_Lagrangian_Force_Data_Insect(path,numSim)

% read in Mag. of Force %
strChoice = 'fX_Lag'; 
fX_Lag =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Tangential Force %
strChoice = 'fY_Lag'; 
fY_Lag = read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Force %
%strChoice = 'fMag'; 
%fLagMag =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Tangential Force %
%strChoice = 'fTan'; 
%fLagNorm = read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);

% read in Mag. of Normal Force %
%strChoice = 'fNorm';
%fLagTan =  read_Force_Scalar_Data_From_vtk(path,numSim,strChoice);
 

clear strChoice first;

