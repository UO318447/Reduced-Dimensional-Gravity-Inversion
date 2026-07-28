%*******************************************************************************
% SYNTHETIC GRAVITY INVERSION - RECTANGULAR PRISMS
%
% Purpose:
%   Runs the synthetic 3D gravity-inversion experiments using RR-GPSO and
%   a rectangular-prism parametrization.
%
% Description:
%   This script generates synthetic gravity data for a prescribed subsurface
%   model composed of rectangular prisms, adds Gaussian noise, configures the
%   RR-GPSO inversion, performs multiple independent inversion runs, and stores
%   the resulting models, predicted gravity responses, convergence information,
%   and execution times.
%
% Main configuration:
%   n_runs      - Number of independent inversion runs.
%   noise_level - Relative noise level added to the synthetic gravity data.
%
% Authors:
%   GPINV Research Group
%   University of Oviedo, Spain
%
% Copyright (c) 2026 GPINV Research Group, University of Oviedo
%
% This repository version preserves the computational workflow used in the
% synthetic gravity-inversion experiments reported in the associated manuscript.
%*******************************************************************************
clearvars
clc
close all
% Reproducibility settings
n_runs = 30;
noise_level = 0.10; % Use 0.05 or 0.10 as reported in the manuscript

for run = 1:n_runs

% -------------------------------------------------------------------------
% 1. SYNTHETIC MODEL AND SEARCH SPACE
% -------------------------------------------------------------------------
% Synthetic bodies: [xc, yc, zc, A, B, C, azimuth, dip, density]
cuerpos_sinteticos = [120 140 -60 35 25 50 30 60 -600; 300 220 -80 60 20 70 90 45 800];

% Lower and upper bounds of the PSO search space
model.lowlimit   = [80  100 -90  40 15 10  0   0 -800  , 240 180 -110 40 10 10  0   0  600];
model.upperlimit = [160 180 -30 100 40 150 180 180 -300, 360 260  -50 90 35 100 180 180 1000];

tic
%==========================================================================
% 2. OBSERVATION GRID
%==========================================================================
pe = 10.0;   % Discretization parameter

% Convert orientation angles from degrees to radians (columns 7 and 8)
cuerpos_sinteticos(:,7:8) = cuerpos_sinteticos(:,7:8)*pi/180;

% Define the observation grid from the spatial extent of the synthetic model
paso = 25.0; 
factor_ancho = 5.0;
ancho = factor_ancho * max(cuerpos_sinteticos(:,4));
limX = [min(cuerpos_sinteticos(:,1))-ancho/2, max(cuerpos_sinteticos(:,1))+ancho/2];
limY = [min(cuerpos_sinteticos(:,2))-ancho/2, max(cuerpos_sinteticos(:,2))+ancho/2];
[x,y] = meshgrid(limX(1):paso:limX(2), limY(1):paso:limY(2));
ptos = [x(:) y(:) zeros(length(x(:)),1)];

rgrid.ptos = ptos;
opfun.ptos = ptos;
opfun.x = x;
opfun.y = y;

%==========================================================================
% 3. SYNTHETIC GRAVITY DATA
%==========================================================================
disp('Generando datos sintéticos: PRISMAS RECTANGULARES')
[gD,~,~,edisc] = GravPrismaRectangular(cuerpos_sinteticos, ptos, 1, 1.0);
gdata = -gD(:,3)*1.0e8;

% Add Gaussian noise to the synthetic gravity data.
data.gobs = gdata + noise_level*randn(size(gdata))*mean(abs(gdata));
%==========================================================================
% 4. INVERSION CONFIGURATION
%==========================================================================
opfun.nparam = 9; % [xc, yc, zc, A, B, C, azimuth, dip, density]
opfun.nstation = length(data.gobs);

% Determine the number of bodies from the model parametrization
ncuerpos = length(model.lowlimit) / opfun.nparam;

disp('Configurando Inversión para: PRISMAS RECTANGULARES')
funobj = @fcost_prisma_rectangular_real;
opfun.npris = ncuerpos; 

%==========================================================================
% 5. RR-GPSO OPTIONS
%==========================================================================
algo_type = 5;
PSO_options; 
options.pso.maxiter = 500; 
options.pso.size = 100; 
options.pso.elitism = 0;
opfun.norm = 2;
opfun.modellog = 0;
opfun.prior.model = []; 
opfun.prior.niter = 2; 
opfun.ordenar = true; 

%==========================================================================
% 6. INVERSE-PROBLEM SOLUTION
%==========================================================================
disp('Iniciando inversión PSO...')
[results] = pso_grav3D(funobj,model,data,options,opfun);

%--------------------------------------------------------------------------
% Extract the best inversion result
%--------------------------------------------------------------------------
[~, jmodel] = min(results.error_hist);
if opfun.modellog == 1
    results.historia = 10.^results.historia;
end
cuerpo_invertido = results.historia(jmodel,:);
if ncuerpos > 1
    cuerpo_invertido = reshape(cuerpo_invertido, [opfun.nparam, ncuerpos])';
end
cuerpo_invertido(:,7:8) = cuerpo_invertido(:,7:8)*pi/180;

%--------------------------------------------------------------------------
% Local refinement of the rectangular-prism model
%--------------------------------------------------------------------------
disp('Refinamiento final (prismas discretizados)...')
mejor_modelo = cuerpo_invertido;
mejor_error = inf;
nref = 5;   
for k = 1:nref
    % Apply a small random perturbation to the current model
    modelo_test = mejor_modelo + 0.02*randn(size(mejor_modelo)).*abs(mejor_modelo);
    % Compute the gravity response of the perturbed discretized model
    [gtest_temp, ~, ~, ~] = GravPrismaRectangular(modelo_test, ptos, 1, 1.0);
    gtest = -gtest_temp(:,3)*1.0e8;
    err = norm(gtest - data.gobs)/norm(data.gobs);
    if err < mejor_error
        mejor_error = err;
        mejor_modelo = modelo_test;
    end
end
cuerpo_invertido = mejor_modelo;

% Compute the final predicted gravity response
[gpre_temp, ~, ~, ediscinv] = GravPrismaRectangular(cuerpo_invertido, ptos, 1, 1.0);
gpre = -gpre_temp(:,3)*1.0e8;
results.gpre = gpre;

Duracion= toc;
% Save the results of the current independent run
save(sprintf('SYN_%s_Run%02d.mat', 'Pris_rect', run), ...
    'results','Duracion','data','gpre', ...
    'cuerpo_invertido','opfun','x','y','ptos');
end
%==========================================================================
% 7. RR-GPSO CONVERGENCE
%==========================================================================
figure('Name', 'PSO convergence')
plot(100*results.error_iter, 'k.-', 'LineWidth', 1.5)
xlabel('Number of iterations')
ylabel('Iteration error (%)')
grid on
% 
%==========================================================================
% 8. ANOMALY MAPS
%==========================================================================
dibuja_anomali(gdata, gpre, x, y)

