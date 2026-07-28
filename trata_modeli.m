function [modelf] = trata_modeli(modeli, opfun)
%*******************************************************************************
% TRATA_MODELI
%
% Purpose:
%   Formats the parameter vector of a candidate model for gravity inversion
%   using multiple geometric bodies and enforces a consistent ordering of
%   the semi-axis parameters.
%
% Syntax:
%   modelf = trata_modeli(modeli, opfun)
%
% Inputs:
%   modeli - Parameter vector representing a candidate model.
%   opfun   - Structure containing the model configuration, including the
%             number of bodies (nelip) and parameters per body (nparam).
%
% Output:
%   modelf  - Formatted model. For multiple bodies, each row contains the
%             parameters of one body. The semi-axis parameters are ordered
%             so that column 4 is greater than or equal to column 5.
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

modelf = modeli;

if opfun.nelip > 1
    modeln = reshape(modeli,[opfun.nparam,opfun.nelip]);
    modelf = modeln';
end

% Enforce a consistent ordering of the semi-axis parameters.
pos = modelf(:,4)<modelf(:,5);
aux = modelf(pos,4);
modelf(pos,4) = modelf(pos,5);
modelf(pos,5) = aux;