%*******************************************************************************
% SYNTHETIC GRAVITY INVERSION EXPERIMENT
%
% Purpose:
%   Runs the synthetic 3D gravity-inversion experiments using RR-GPSO and
%   reduced-dimensional geometric parametrizations.
%
% Description:
%   This script generates synthetic gravity data for a prescribed subsurface
%   model, adds Gaussian noise, configures the RR-GPSO inversion, performs
%   multiple independent inversion runs, and stores the resulting models,
%   predicted gravity responses, convergence information, and execution times.
%
%   The script supports square-based prism and prolate-ellipsoid
%   parametrizations through the variables tipo_sintetico and tipo_inversion.
%
% Main configuration:
%   n_runs         - Number of independent inversion runs.
%   noise_level    - Relative noise level added to the synthetic gravity data.
%   tipo_sintetico - Geometry used to generate the synthetic data.
%   tipo_inversion - Geometry used for the inversion.
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
    %--------------------------------------------------------------------------
    % 1. GEOMETRIC PARAMETRIZATION
    %--------------------------------------------------------------------------
    % 1 = Square-based prism | 2 = Prolate ellipsoid
    tipo_sintetico = 2; % Geometry used to generate the synthetic data
    tipo_inversion = 2; % Geometry used for the inversion

    %--------------------------------------------------------------------------
    % 2. SYNTHETIC MODEL AND SEARCH SPACE
    %--------------------------------------------------------------------------
    % Synthetic bodies: [xc, yc, zc, A, B, azimuth, dip, density]
    cuerpos_sinteticos = [120 140 -60 35 25 30 60 -600; 300 220 -80 60 20 90 45 800];

    % Lower and upper bounds of the PSO search space
    model.lowlimit   = [80  100 -90  40 15   0   0 -800, 240 180 -110 40 10   0   0  600];
    model.upperlimit = [160 180 -30 100 40 180 180 -300, 360 260  -50 90 35 180 180 1000];


    tic
    %==========================================================================
    % 3. OBSERVATION GRID
    %==========================================================================
    pe = 10.0;   % Discretization parameter

    % Convert orientation angles from degrees to radians
    cuerpos_sinteticos(:,6:7) = cuerpos_sinteticos(:,6:7)*pi/180;

    % Define the observation grid from the spatial extent of the synthetic model
    paso = 25.0; % Observation-grid spacing
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
    % 4. SYNTHETIC GRAVITY DATA
    %==========================================================================
    if tipo_sintetico == 1
        disp('Generando datos sintéticos: PRISMAS')
        [gD,~,~,edisc] = GravPrismaCuadrado(cuerpos_sinteticos, ptos, 1, 1.0);
        gdata = -gD(:,3)*1.0e8;
    elseif tipo_sintetico == 2
        disp('Generando datos sintéticos: ELIPSOIDES PROLATOS')
        [opfun_temp, gdata] = problema_directo_prolato(rgrid, cuerpos_sinteticos);
        edisc = opfun_temp.edisc;
    end

    % Add Gaussian noise to the synthetic gravity data.
    data.gobs = gdata + noise_level*randn(size(gdata))*mean(abs(gdata));

    %==========================================================================
    % 5. INVERSION CONFIGURATION
    %==========================================================================
    opfun.nparam = 8;
    opfun.nstation = length(data.gobs);

    % Determine the number of geometric bodies from the model parametrization.
    ncuerpos = length(model.lowlimit) / opfun.nparam;

    if tipo_inversion == 1
        disp('Configurando Inversión para: PRISMAS CUADRADOS')
        funobj = @fcost_prisma_cuadrado_real;
        opfun.npris = ncuerpos; % Number of prisms
    elseif tipo_inversion == 2
        disp('Configurando Inversión para: ELIPSOIDES PROLATOS')
        funobj = @fcost_elipsoide_prolato_real;
        opfun.nelip = ncuerpos; % Number of ellipsoids
    end

    %==========================================================================
    % 6. RR-GPSO OPTIONS
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
    % 7. INVERSE-PROBLEM SOLUTION
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
    cuerpo_invertido(:,6:7) = cuerpo_invertido(:,6:7)*pi/180;

    % Compute the predicted gravity response for the selected parametrization

    if tipo_inversion == 1

        %--------------------------------------------------------------------------
        % Local refinement for the square-based prism parametrization
        %--------------------------------------------------------------------------
        disp('Refinamiento final (prismas discretizados)...')

        mejor_modelo = cuerpo_invertido;
        mejor_error = inf;

        nref = 5;   % Number of local refinement trials

        for k = 1:nref

            % Apply a small random perturbation to the current model
            modelo_test = mejor_modelo + 0.02*randn(size(mejor_modelo)).*abs(mejor_modelo);

            % Compute the gravity response of the perturbed discretized model
            [gtest_temp, ~, ~, ~] = GravPrismaCuadrado(modelo_test, ptos, 1, 1.0);
            gtest = -gtest_temp(:,3)*1.0e8;

            err = norm(gtest - data.gobs)/norm(data.gobs);

            if err < mejor_error
                mejor_error = err;
                mejor_modelo = modelo_test;
            end
        end

        cuerpo_invertido = mejor_modelo;

        % Compute the final predicted gravity response
        [gpre_temp, ~, ~, ediscinv] = GravPrismaCuadrado(cuerpo_invertido, ptos, 1, 1.0);

        gpre = -gpre_temp(:,3)*1.0e8;

    elseif tipo_inversion == 2
        % Enforce a consistent ordering of the prolate-ellipsoid semi-axes
        pos = cuerpo_invertido(:,4) < cuerpo_invertido(:,5);
        aux = cuerpo_invertido(pos,4);
        cuerpo_invertido(pos,4) = cuerpo_invertido(pos,5);
        cuerpo_invertido(pos,5) = aux;

        [gpre_temp, ediscinv] = GravElipProlato(cuerpo_invertido, ptos);
        gpre = -gpre_temp(:,3)*1.0e8;
    end

    results.gpre = gpre;
    Duracion=toc

    % Save the results of the current independent run.
    if tipo_inversion == 1
        pD = 'SqPrism';
    elseif tipo_inversion == 2
        pD = 'Elip';
    end
    save(sprintf('SYN_%s_Run%02d.mat', pD, run), ...
        'results','Duracion','data','gpre', ...
        'cuerpo_invertido','opfun','x','y','ptos');
end

%==========================================================================
% 8. RESULTS VISUALIZATION
%==========================================================================
disp('Generando gráficas...')

%--------------------------------------------------------------------------
% 8.1. Generate the geometries required for 3D visualization
%--------------------------------------------------------------------------

% Synthetic bodies
if tipo_sintetico == 1
    [~, ~, ~, edisc] = GravPrismaCuadrado(cuerpos_sinteticos, ptos);
elseif tipo_sintetico == 2
    % Obtain the surface representation of the synthetic ellipsoids.
    [opfun_temp, ~] = problema_directo_prolato(rgrid, cuerpos_sinteticos);
    edisc = opfun_temp.edisc;
end

% Inverted bodies
if tipo_inversion == 1
    [~, ~, ~, ediscinv] = GravPrismaCuadrado(cuerpo_invertido, ptos);
elseif tipo_inversion == 2
    [~, ediscinv] = GravElipProlato(cuerpo_invertido, ptos);
end

%--------------------------------------------------------------------------
% 8.2. 3D comparison of the synthetic and inverted models
%--------------------------------------------------------------------------
figure
hold on

h_sint = [];
for i = 1:length(edisc)
    if isempty(edisc{i}) || size(edisc{i},2) ~= 3
        continue
    end
    C_sintetico = convhull(edisc{i}(:,1), edisc{i}(:,2), edisc{i}(:,3));
    h_sint = trisurf(C_sintetico, ...
        edisc{i}(:,1), edisc{i}(:,2), edisc{i}(:,3), ...
        'FaceColor', 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
end

h_inv = [];

for j = 1:length(ediscinv)
    if isempty(ediscinv{j}) || size(ediscinv{j},2) ~= 3
        continue
    end
    C_invertido = convhull(ediscinv{j}(:,1), ediscinv{j}(:,2), ediscinv{j}(:,3));
    h_inv = trisurf(C_invertido, ...
        ediscinv{j}(:,1), ediscinv{j}(:,2), ediscinv{j}(:,3), ...
        'FaceColor', 'r', 'FaceAlpha', 0.8, 'EdgeColor', 'k');
end

% Observation points at the surface
scatter3(ptos(:,1), ptos(:,2), ptos(:,3), 5, 'k', 'filled');

% Configure the 3D visualization
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Profundidad Z (m)')
legend([h_sint, h_inv], {'Synthetic Model', 'Inverted Model'}, 'Location', 'best')
grid on; view(3); axis equal;
hold off

%--------------------------------------------------------------------------
% 8.3. RR-GPSO convergence
%--------------------------------------------------------------------------
figure('Name', 'PSO convergence')
plot(100*results.error_iter, 'k.-', 'LineWidth', 1.5)
xlabel('Number of iterations')
ylabel('Iteration error (%)')
grid on

%--------------------------------------------------------------------------
% 8.4. Observed and predicted gravity-anomaly maps
%--------------------------------------------------------------------------
dibuja_anomali(gdata, gpre, x, y)



