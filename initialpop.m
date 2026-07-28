function [parent]=initialpop(talla,lowlimit,upperlimit);
%*******************************************************************************
% INITIALPOP
%
% Purpose:
%   Generates a random initial population within the prescribed search-space
%   bounds.
%
% Syntax:
%   parent = initialpop(talla, lowlimit, upperlimit)
%
% Inputs:
%   talla      - Number of models (particles) in the population.
%   lowlimit   - Lower bounds of the model parameters.
%   upperlimit - Upper bounds of the model parameters.
%
% Output:
%   parent     - Randomly generated population. Each row represents one model
%                within the interval defined by lowlimit and upperlimit.
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
% Random initial population in [lowlimit, upperlimit]
nparam=length(lowlimit);
rango=upperlimit-lowlimit;
parent=(ones(talla,1)*rango).*(rand(talla,nparam))+(ones(talla,1)*lowlimit);