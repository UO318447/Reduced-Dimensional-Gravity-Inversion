function [modelf] = trata_modelip(modeli,opfun)
%*******************************************************************************
% TRATA_MODELIP
%
% Purpose:
%   Formats the parameter vector of a candidate prism model according to
%   the number of prisms and, when explicitly requested, orders the prisms
%   by their first coordinate.
%
% Syntax:
%   modelf = trata_modelip(modeli, opfun)
%
% Inputs:
%   modeli - Parameter vector representing a candidate prism model.
%   opfun   - Structure containing the model configuration, including the
%             number of prisms (npris), the number of parameters per prism
%             (nparam), and the optional ordering flag (ordenar).
%
% Output:
%   modelf  - Formatted model in which each row contains the parameters
%             of one prism.
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

if opfun.npris > 1
    modeln = reshape(modeli,[opfun.nparam,opfun.npris]);
    modelf = modeln';
else
    modelf = reshape(modeli,[opfun.nparam,1])';
end

% Sort the prisms only when explicitly requested
if opfun.npris > 1 && isfield(opfun,'ordenar') && opfun.ordenar
    [~,idx] = sort(modelf(:,1));
    modelf = modelf(idx,:);
end
