function [misfit,swarm] = fcost_prisma_rectangular_real(swarm,data,opciones,opfun)
%*******************************************************************************
% FCOST_PRISMA_RECTANGULAR_REAL
%
% Purpose:
%   Evaluates the normalized gravity-data misfit for models parametrized
%   as rectangular prisms.
%
% Syntax:
%   [misfit, swarm] = fcost_prisma_rectangular_real(swarm, data, opciones, opfun)
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
%   trata_modelip.m
%   GravPrismaRectangular.m
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

for i = 1:size(swarm,1)

    modeli = swarm(i,:);
    pris = trata_modelip(modeli,opfun);
    
    % Reorder the model vector for PSO.
    B = pris';
    Bf = B(:);
    swarm(i,:) = Bf';
    
    %-----------------------------------------
    % Convert orientation angles from degrees to radians
    %-----------------------------------------
    pris(:,7:8) = pris(:,7:8)*pi/180;
    
    %-----------------------------------------
    % Overlap constraint for rectangular prisms
    %-----------------------------------------
    invalid = false;

    if opfun.npris > 1
        for k = 1:opfun.npris-1
            for j = k+1:opfun.npris
                
                dx = abs(pris(k,1)-pris(j,1));
                dy = abs(pris(k,2)-pris(j,2));
                dz = abs(pris(k,3)-pris(j,3));
                
                % Prism half-dimensions
                ax = pris(k,4);
                ay = pris(k,5);
                az = pris(k,6);
                
                bx = pris(j,4);
                by = pris(j,5);
                bz = pris(j,6);
                
                overlap_x = dx < (ax + bx);
                overlap_y = dy < (ay + by);
                overlap_z = dz < (az + bz);
                
                if overlap_x && overlap_y && overlap_z
                    invalid = true;
                    break
                end
            end
            if invalid
                break
            end
        end
    end
    
    % Penalize candidate models containing overlapping prisms
    if invalid
        misfit(i) = 1e6;
        continue
    end
    
    %-----------------------------------------
    % Forward gravity calculation
    %-----------------------------------------
    gD = GravPrismaRectangular(pris, ptos, 1, 1.0);
    
    % Extract the vertical gravity component
    gcal = -gD(:,3) * 1.0e8;
    
    %-----------------------------------------
    % Normalized data misfit
    %-----------------------------------------
    misfit(i) = norm(gcal(:)-data.gobs(:),opfun.norm) / ...
                norm(data.gobs(:),opfun.norm);

end
end