function [good_models] = posterior(opfun,options,results,Etol)
%*******************************************************************************
% POSTERIOR MODEL SELECTION
%
% Purpose:
%   Selects representative models from the RR-GPSO search history according
%   to the swarm dispersion and a prescribed error tolerance.
%
% Inputs:
%   opfun   - Problem-dependent parameter structure.
%   options - RR-GPSO configuration structure.
%   results - Structure containing the inversion results and search history.
%   Etol    - Error tolerance used to select acceptable models.
%
% Output:
%   good_models - Selected models for posterior analysis.
%
% Authors:
%   GPINV Research Group
%   University of Oviedo, Spain
%
% Copyright (c) 2026 GPINV Research Group, University of Oviedo
%*******************************************************************************

Edis = 10;

disper = results.stat.distmed;
disper = disper/disper(1)*100;

% Identify iterations with relative swarm dispersion above the threshold
iterd = find(disper >= Edis);

nparticles = options.pso.size;
niter = options.pso.maxiter;
good_models = [];

% Select models satisfying the error tolerance when swarm dispersion
% remains above the prescribed threshold
for i = 1:length(iterd)
    inicio = 1 + (iterd(i)-1)*nparticles;
    fin = inicio + nparticles - 1;
    igoodm = find(results.error_hist(inicio:fin) <= Etol);
    good_models = [good_models; results.historia(igoodm,:)];
end

% For iterations below the dispersion threshold, retain the median model
% as a representative solution
rest_iter = setdiff(1:niter,iterd);

for i = 1:length(rest_iter)
    inicio = 1 + (rest_iter(i)-1)*nparticles;
    fin = inicio + nparticles - 1;
    goodm = median(results.historia(inicio:fin,:));
    good_models = [good_models; goodm];
end