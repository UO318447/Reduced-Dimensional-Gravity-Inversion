for i = 1:30 % Acaba en la última linea
%*******************************************************************************
% REAL-DATA 3D GRAVITY INVERSION USING RR-GPSO
%
% Purpose:
%   Performs 3D gravity inversion of real gravity data using RR-GPSO and
%   alternative reduced-dimensional geometric parametrizations.
%
% Description:
%   This script loads a real gravity dataset, defines the corresponding
%   search space, configures the selected geometric parametrization, performs
%   multiple independent RR-GPSO inversion runs, and stores the resulting
%   models, predicted gravity responses, convergence information, and
%   execution times.
%
% Supported parametrizations:
%   1 - Square-based prisms
%   2 - Prolate ellipsoids
%   3 - Rectangular prisms
%
% Authors:
%   GPINV Research Group
%   University of Oviedo, Spain
%
% Copyright (c) 2026 GPINV Research Group, University of Oviedo
%
% This repository version preserves the computational workflow used in the
% real-data gravity-inversion experiments reported in the associated manuscript.
%*******************************************************************************
clear all
clc
close all
% -------------------------------------------------------------------------
% 1. CONTROL PANEL: DATASET AND GEOMETRIC PARAMETRIZATION
% -------------------------------------------------------------------------
% Select the dataset: 2 = Nirano | 3 = California 
dataset_choice = 3; 

% Select the inversion geometry:
% 1 = Square-based prisms | 2 = Prolate ellipsoids | 3 = Rectangular prisms
inversion_choice = 3;

% -------------------------------------------------------------------------
% 2. GLOBAL RR-GPSO OPTIONS
% -------------------------------------------------------------------------
algo_type = 5;                 
options.pso.maxiter = 300;     % Maximum number of iterations
options.pso.size = 150;        % Maximum number of iterations
options.pso.elitism = 0;       % Elitism option

%==========================================================================
% 3. DATA LOADING AND SEARCH-SPACE DEFINITION
%==========================================================================
if dataset_choice == 1
    error('Dataset 1 is not included in this repository version.');
elseif dataset_choice == 2
    disp('>>> Cargando datos de NIRANO (VOLCÁN)...')
    datos = xlsread('Volcan-1.xlsx'); 
    x = datos(:,1); y = datos(:,2); z = datos(:,3); 
    gobs = datos(:,4); % Gravity observations in microGal

% Search space for the four-body Nirano model.
    model.lowlimit   = [min(x) min(y) -250 10   10  0  0 -500, min(x) min(y) -250 10   10  0  0 -500, min(x) min(y) -250 10   10  0  0  0, min(x) min(y) -250 10  10  0 0  0]; 
    model.upperlimit = [max(x) max(y)   0 500 300 180 180  0,    max(x) max(y)    0 500 300 180 180 0,    max(x) max(y)   0 500 300 90 180 500, max(x) max(y) 0 500 300 180 180 500];

elseif dataset_choice == 3
    disp('>>> Cargando datos de CALIFORNIA (Clear Lake Volcanic Field)...')
    datos = xlsread('CaliforniaISO.xlsx');
    x = datos(:,1); y = datos(:,2); z = datos(:,3); 
    gobs = datos(:,4); 
    
    % Reduce the spatial density of observation stations.
dmin = 3000; % Minimum distance between retained stations (m)
    idx = 1;
    for i = 2:length(x)
        if sqrt((x(i)-x(idx(end)))^2 + (y(i)-y(idx(end)))^2) > dmin
            idx(end+1) = i;
        end
    end
    x = x(idx); y = y(idx); z = z(idx); gobs = gobs(idx);
    
    disp(['    -> Estaciones originales: ', num2str(length(datos(:,1)))]);
    disp(['    -> Estaciones reducidas: ', num2str(length(gobs))]);
    
    % RR-GPSO configuration for the California dataset
    options.pso.maxiter = 300; 
    options.pso.size = 100;     
    
    
%==========================================================================
% SEARCH SPACE FOR THE 12-BODY CALIFORNIA MODEL
%==========================================================================

low1 = [513005-2000 4282550-2000 -2500   3000   1500   0   0  -600];
up1  = [513005+2000 4282550+2000  -200    9000   3000 180  20    0];

low2 = [519445-2000 4341459-2000 -8000 1000 800   0  0 -600];
up2  = [519445+2000 4341459+2000 -500  4000 3000 180 20 0];

low3 = [520641-2000 4262962-2000 -4500  100  500  0  0 -600];
up3  = [520641+2000 4262962+2000  -500  8000 2500 180 20 0];

low4 = [505259-2000 4337645-2000 -5000 2000  1000   0   0 -600];
up4  = [505259+2000 4337645+2000 -1000 8000  3000 180 20 0];

low5 = [535800-2000 4272473-2000 -7000 2000 800 0 0 -600];
up5  = [535800+2000 4272473+2000 -1000 8000 3000 180 20 0];

low6 = [513473-2000 4292731-2000 -8000 2000 1000 0 0 -600];
up6  = [513473+2000 4292731+2000 -2000  8000 4000 180 20 -50];

low7 = [538684-2000 4328919-2000   -9000    1000   1000   0   0  -800];
up7  = [538684+2000 4328919+2000   -2000   20000   5000 180  20  -100];

low8 = [520794-2000 4305304-2000  -8000  10000  1000   0    0  -800];
up8  = [520794+2000 4305304+2000  -1000  20000  4000 180  20   800];

low9 = [534585-2000 4253196-2000  -5000   100   100    0    0   -1000];
up9  = [534585+2000 4253196+2000   -200  3000  3000  180   20      0];

low10 = [540847-2000 4315226-2000   -9000   100   100   0   0  -300];
up10  = [540847+2000 4315226+2000   -2000  3000 3000 180  20  -100];

low11 = [min(x) min(y)  -8000  1000  1000   0    0  -800];
up11  = [max(x) max(y)  -1000  20000  4000 180  180   800];

low12 = [min(x) min(y)  -8000  1000  1000   0    0  -800];
up12  = [max(x) max(y)  -1000  20000  4000 180  180   800];

model.lowlimit   = [low1, low2, low3, low4, low5, low6, low7, low8, low9, low10, low11, low12];
model.upperlimit = [up1,  up2,  up3,  up4,  up5,  up6,   up7,  up8,  up9, up10, up11, up12];
    
end

%--------------------------------------------------------------------------
% Automatic adaptation of the search space for rectangular prisms
%--------------------------------------------------------------------------
% If the bounds are defined using eight parameters per body, an additional
% semi-axis parameter C is introduced using the bounds defined for B
if inversion_choice == 3 && isfield(model, 'lowlimit') && mod(length(model.lowlimit), 8) == 0
    nc = length(model.lowlimit) / 8;
    L_mat = reshape(model.lowlimit, [8, nc]);
    U_mat = reshape(model.upperlimit, [8, nc]);
    
    % Introduce parameter C using the bounds of parameter B
    L_new = [L_mat(1:5, :); L_mat(5, :); L_mat(6:8, :)];
    U_new = [U_mat(1:5, :); U_mat(5, :); U_mat(6:8, :)];
    
    model.lowlimit = L_new(:)';
    model.upperlimit = U_new(:)';
    disp('>>> (Auto-adaptación): Límites de parámetro C añadidos automáticamente clonando los de B.');
end

% Store observation coordinates and gravity data
ptos = [x y z];
data.gobs = gobs;
opfun.ptos = ptos; opfun.x = x; opfun.y = y; opfun.z = z;

%==========================================================================
% 4. INVERSE-MODEL CONFIGURATION
%==========================================================================
opfun.nstation = length(data.gobs);

% Set the number of parameters per body according to the geometry.
if inversion_choice == 3
    opfun.nparam = 9;
else
    opfun.nparam = 8;
end

ncuerpos = length(model.lowlimit) / opfun.nparam; % Number of model bodies

if inversion_choice == 1
    disp(['>>> Configurando inversión: ', num2str(ncuerpos), ' PRISMAS CUADRADOS'])
    funobj = @fcost_prisma_cuadrado_real;
    opfun.npris = ncuerpos;
elseif inversion_choice == 2
    disp(['>>> Configurando inversión: ', num2str(ncuerpos), ' ELIPSOIDES PROLATOS'])
    funobj = @fcost_elipsoide_prolato_real;
    opfun.nelip = ncuerpos;
elseif inversion_choice == 3
    disp(['>>> Configurando inversión: ', num2str(ncuerpos), ' PRISMAS RECTANGULARES'])
    funobj = @fcost_prisma_rectangular_real;
    opfun.npris = ncuerpos; 
end

% Internal options required by the objective function.
PSO_options; 
opfun.norm = 1; opfun.modellog = 0; opfun.prior.model = []; opfun.prior.niter = 2; 

%==========================================================================
% 5. INVERSE-PROBLEM SOLUTION
%==========================================================================
disp('>>> Iniciando optimización PSO...')
tic
[results] = pso_grav3D(funobj, model, data, options, opfun);
Duracion=toc

% Extract the best inversion result
[~, jmodel] = min(results.error_hist);
if opfun.modellog == 1
    results.historia = 10.^results.historia;
end
cuerpo_invertido = results.historia(jmodel,:);
if ncuerpos > 1
    cuerpo_invertido = reshape(cuerpo_invertido, [opfun.nparam, ncuerpos])';
end

% Convert orientation angles from degrees to radians.
if inversion_choice == 3
    cuerpo_invertido(:,7:8) = cuerpo_invertido(:,7:8)*pi/180; 
else
    cuerpo_invertido(:,6:7) = cuerpo_invertido(:,6:7)*pi/180; 
end

% Compute the predicted gravity response and 3D geometry
if inversion_choice == 1
    [gpre_temp, ~, ~, ediscinv] = GravPrismaCuadrado(cuerpo_invertido, ptos, 1, 1.0, 0.0);
elseif inversion_choice == 2
    pos = cuerpo_invertido(:,4) < cuerpo_invertido(:,5); % Enforce a consistent ordering of the prolate-ellipsoid semi-axes
    aux = cuerpo_invertido(pos,4);
    cuerpo_invertido(pos,4) = cuerpo_invertido(pos,5);
    cuerpo_invertido(pos,5) = aux;
    
    [gpre_temp, ediscinv] = GravElipProlato(cuerpo_invertido, ptos);
elseif inversion_choice == 3
    [gpre_temp, ~, ~, ediscinv] = GravPrismaRectangular(cuerpo_invertido, ptos, 1, 1.0);
end

gpre = -gpre_temp(:,3)*1.0e8; % Convert from m/s^2 to microGal
results.gpre = gpre;

% Generate the output filename
if dataset_choice == 1; prefijo_datos = 'Cerr';
elseif dataset_choice == 2; prefijo_datos = 'Nirano';
elseif dataset_choice == 3; prefijo_datos = 'Cali'; end

if inversion_choice == 1; prefijo_geom = 'Pris';
elseif inversion_choice == 2; prefijo_geom = 'Elip';
elseif inversion_choice == 3; prefijo_geom = 'PRect'; end

error_final = 100 * min(results.error_hist); 
str_error = sprintf('Err%.2f', error_final); 
str_error = strrep(str_error, '.', 'p');    
timestamp = num2str(Duracion); 
timestamp = strrep(timestamp, '.', 'p');
nombre_archivo = [prefijo_datos, '_', prefijo_geom, '_', str_error, '_', timestamp];

disp(['>>> Guardando resultados como: ', nombre_archivo, '.mat']);
save(nombre_archivo);

%%
%==========================================================================
% 6. POSTERIOR ANALYSIS
%==========================================================================
disp('>>> Calculando análisis a posteriori...')
Etol = 0.20; 
[good_models] = posterior(opfun, options, results, Etol);

%==========================================================================
% 7. RESULTS VISUALIZATION
%==========================================================================
disp('>>> Generando gráficas universales...')

% 7.1. Inverted 3D model
f1 = figure('Name', 'Inverted 3D Model', 'Position', [100, 100, 800, 600]);
hold on
nelinv = length(ediscinv);
colores = lines(nelinv); 
h_cuerpos = []; 
nombres_cuerpos = {}; 
for i=1:nelinv
    Cinv = convhull(ediscinv{i}(:,1), ediscinv{i}(:,2), ediscinv{i}(:,3));
    h_cuerpos(i) = trisurf(Cinv, ediscinv{i}(:,1), ediscinv{i}(:,2), ediscinv{i}(:,3), ...
        'FaceColor', colores(i,:), ... 
        'FaceAlpha', 0.4, ...          
        'EdgeColor', 'none');          
    nombres_cuerpos{i} = ['Body ', num2str(i)]; 
end
h_estaciones = scatter3(x, y, z, 15, 'k', 'filled'); 
title(['3D Model: ', num2str(ncuerpos), ' Bodies']); 
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)')
grid on; view(3); axis equal; 
legend([h_cuerpos, h_estaciones], [nombres_cuerpos, {'Stations (Topography)'}], 'Location', 'northeast');

ax = gca;
ax.Position = [0.05 0.1 0.62 0.8]; 
texto_datos = {'INVERTED MODEL', '-------------------'}; 
for i = 1:ncuerpos
    texto_datos{end+1} = sprintf('Body %d:', i);
    texto_datos{end+1} = sprintf(' Center : [%.0f, %.0f, %.0f]', cuerpo_invertido(i,1), cuerpo_invertido(i,2), cuerpo_invertido(i,3));

    if inversion_choice == 3
        texto_datos{end+1} = sprintf(' Axes   : [%.0f, %.0f, %.0f]', cuerpo_invertido(i,4), cuerpo_invertido(i,5), cuerpo_invertido(i,6));
        texto_datos{end+1} = sprintf(' Angles: [%.1f°, %.1f°]', cuerpo_invertido(i,7)*180/pi, cuerpo_invertido(i,8)*180/pi);
        texto_datos{end+1} = sprintf(' Density: %.0f kg/m^3', cuerpo_invertido(i,9));
    else
        texto_datos{end+1} = sprintf(' Axes   : [%.0f, %.0f]', cuerpo_invertido(i,4), cuerpo_invertido(i,5));
        texto_datos{end+1} = sprintf(' Angles: [%.1f°, %.1f°]', cuerpo_invertido(i,6)*180/pi, cuerpo_invertido(i,7)*180/pi);
        texto_datos{end+1} = sprintf(' Density: %.0f kg/m^3', cuerpo_invertido(i,8));
    end
    texto_datos{end+1} = ' '; 
end
annotation('textbox', [0.68 0.1 0.30 0.8], 'String', texto_datos, 'EdgeColor', 'k', ...
    'BackgroundColor', [1 1 1], 'FaceAlpha', 0.9, 'FitBoxToText', 'on', 'FontSize', 8, 'FontName', 'Courier');
hold off

% 7.2. Observed and predicted gravity-anomaly maps
num_pts = 100; 
[X_grid, Y_grid] = meshgrid(linspace(min(x), max(x), num_pts), linspace(min(y), max(y), num_pts));
G_obs_grid = griddata(x, y, data.gobs, X_grid, Y_grid, 'natural');
G_pre_grid = griddata(x, y, gpre, X_grid, Y_grid, 'natural');

f2 = figure('Name', 'Anomaly map', 'Position', [150, 150, 1000, 400]);
subplot(1, 2, 1)
pcolor(X_grid, Y_grid, G_obs_grid); shading interp; hold on;
scatter(x, y, 5, 'k', 'filled'); hold off;
title('Observed anomaly (\muGal)'); xlabel('X (m)'); ylabel('Y (m)');
colorbar; colormap(jet); axis equal tight;

subplot(1, 2, 2)
pcolor(X_grid, Y_grid, G_pre_grid); shading interp; hold on;
scatter(x, y, 5, 'k', 'filled'); hold off;
title('Predicted anomaly (\muGal)'); xlabel('X (m)'); ylabel('Y (m)');
colorbar; colormap(jet); axis equal tight;

% Use consistent color limits for observed and predicted anomalies
if dataset_choice == 2 % Nirano
    vmin = -1500; vmax = 1500;
elseif dataset_choice == 3 % California
    vmin = min(data.gobs); vmax = max(data.gobs);
end
subplot(1,2,1); clim([vmin vmax]); 
subplot(1,2,2); clim([vmin vmax]);

% 7.3. Observed and predicted gravity responses
f3 = figure('Name', 'Ajuste y Convergencia', 'Position', [200, 200, 1000, 400]);
plot(data.gobs, 'k.-', 'DisplayName', 'Observed', 'LineWidth', 1.5)
hold on; plot(gpre, 'r.-', 'DisplayName', 'Predicted', 'LineWidth', 1.5); hold off;
title('Ajuste por Estación'); xlabel('Station number'); ylabel('Gravity (\muGal)');
legend('Location', 'northwest'); grid on;

% ==========================================================================
% 8. SAVE GENERATED FIGURES
% ==========================================================================
nombre_carpeta = [prefijo_datos, '_', prefijo_geom, '_', str_error];
nombre_base = [prefijo_datos, '_', prefijo_geom,'_', str_error];

if ~exist(nombre_carpeta, 'dir')
    mkdir(nombre_carpeta);
end

ruta_f1 = fullfile(nombre_carpeta, [nombre_base, '_Modelo3D']);
saveas(f1, [ruta_f1, '.png']); 
saveas(f1, [ruta_f1, '.fig']); 

ruta_f2 = fullfile(nombre_carpeta, [nombre_base, '_MapasAnomalia']);
saveas(f2, [ruta_f2, '.png']);
saveas(f2, [ruta_f2, '.fig']);

ruta_f3 = fullfile(nombre_carpeta, [nombre_base, '_Convergencia']);
saveas(f3, [ruta_f3, '.png']);
saveas(f3, [ruta_f3, '.fig']);
end