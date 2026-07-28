function [misfit, swarm] = fcost_elipsoide_prolato_real(swarm, data, opciones, opfun)
%*******************************************************************************
% FCOST_ELIPSOIDE_PROLATO_REAL
%
% Purpose:
%   Evaluates the normalized gravity-data misfit for models parametrized
%   as prolate ellipsoids.
%
% Syntax:
%   [misfit, swarm] = fcost_elipsoide_prolato_real(swarm, data, opciones, opfun)
%
% Inputs:
%   swarm     - Population of candidate models. Each row represents one model.
%   data      - Structure containing the observed gravity data (data.gobs).
%   opciones  - Structure containing inversion options.
%   opfun     - Structure containing problem-dependent parameters, including
%               observation points, norm definition, and model configuration.
%
% Outputs:
%   misfit    - Normalized gravity-data misfit for each candidate model.
%   swarm     - Population after internal model formatting.
%
% Dependencies:
%   trata_modeli.m
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
ptos = opfun.ptos;
for i=1:size(swarm,1)
    modeli = swarm(i,:);
    [elip2] = trata_modeli(modeli,opfun);
    
    B = elip2';
    Bf = B(:);
    swarm(i,:) = Bf';
    
    elip2(:,6:7) = elip2(:,6:7)*pi/180;
    
    gD = GravElipProlato(elip2,ptos);
    
    % Sum the vertical gravity contribution of all prolate ellipsoids.
    gcal = -sum(gD(:, 3:3:end), 2) * 1.0e8;
     
    misfit(i) = norm(gcal(:)-data.gobs(:),opfun.norm)/norm(data.gobs(:),opfun.norm);
    rms(i) = norm(gcal(:)-data.gobs(:),opfun.norm)/sqrt(length(data.gobs));
end