function [opfun,data] = problema_directo_prolato(rgrid,elip)
%*******************************************************************************
% PROBLEMA_DIRECTO_PROLATO
%
% Purpose:
%   Computes the forward gravity response of a model parametrized by
%   prolate ellipsoids at the specified observation points.
%
% Syntax:
%   [opfun, data] = problema_directo_prolato(rgrid, elip)
%
% Inputs:
%   rgrid - Structure containing the observation-point coordinates (ptos).
%   elip  - Matrix defining the prolate ellipsoid model.
%
% Outputs:
%   opfun - Structure containing auxiliary information returned by the
%           forward-modeling routine.
%   data  - Vertical gravity response expressed in microGal.
%
% Dependencies:
%   GravElipProlato.m
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

[grave,edisc] = GravElipProlato(elip,rgrid.ptos);
opfun.edisc = edisc;

% Convert the vertical gravity component to microGal and adjust its sign
data = -grave(:,3)*1.0e8;