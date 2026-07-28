function [vox,dvol] = DiscretizaPrisma(a,b,pe)
%*******************************************************************************
% DISCRETIZAPRISMA
%
% Purpose:
%   Discretizes a square-based rectangular prism into a regular grid of
%   volumetric elements for numerical gravity calculations.
%
% Syntax:
%   [vox, dvol] = DiscretizaPrisma(a, b, pe)
%
% Inputs:
%   a    - Half-length of the square base in the X and Y directions.
%   b    - Half-length of the prism in the Z direction.
%   pe   - Parameter controlling the spatial discretization resolution.
%
% Outputs:
%   vox  - Matrix containing the coordinates of the voxel centers.
%   dvol - Volume of each voxel.
%
% Authors:
%   GPINV Research Group
%   University of Oviedo, Spain
%
% Copyright (c) 2026 GPINV Research Group, University of Oviedo
%
% This repository version preserves the computational logic used in the
% associated gravity-inversion experiments.
%*******************************************************************************

% Number of divisions used for the regular discretization.
n = max(4, round(6/pe));

dx = 2*a / n;
dy = 2*a / n;
dz = 2*b / n;

% Coordinates of the voxel centers.
x = linspace(-a + dx/2, a - dx/2, n);
y = linspace(-a + dy/2, a - dy/2, n);
z = linspace(-b + dz/2, b - dz/2, n);

[X,Y,Z] = meshgrid(x,y,z);

vox = [X(:), Y(:), Z(:)];

% Volume of each voxel.
dvol = dx * dy * dz;

end