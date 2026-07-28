function [vox,dvol] = DiscretizaPrismaRectangular(a,b,c,pe)
%*******************************************************************************
% DISCRETIZAPRISMARECTANGULAR
%
% Purpose:
%   Discretizes a rectangular prism into a regular grid of volumetric
%   elements for numerical gravity calculations.
%
% Syntax:
%   [vox, dvol] = DiscretizaPrismaRectangular(a, b, c, pe)
%
% Inputs:
%   a    - Half-dimension of the prism in the X direction.
%   b    - Half-dimension of the prism in the Y direction.
%   c    - Half-dimension of the prism in the Z direction.
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

n = max(4, round(6/pe));
Ntarget = n^3;

ratio = [a b c] / max([a b c]);

nx = round(n * ratio(1));
ny = round(n * ratio(2));
nz = round(n * ratio(3));

% Scale the discretization to control the total number of voxels
scale = (nx * ny * nz / Ntarget)^(1/3);

if scale > 1
    nx = round(nx / scale);
    ny = round(ny / scale);
    nz = round(nz / scale);
end

% Enforce a minimum number of divisions along each dimension
nx = max(nx,4); 
ny = max(ny,4); 
nz = max(nz,4);

dx = 2*a / nx;
dy = 2*b / ny;
dz = 2*c / nz;

x = linspace(-a + dx/2, a - dx/2, nx);
y = linspace(-b + dy/2, b - dy/2, ny);
z = linspace(-c + dz/2, c - dz/2, nz);

[X,Y,Z] = meshgrid(x,y,z);
vox = [X(:), Y(:), Z(:)];

dvol = dx * dy * dz;

end



