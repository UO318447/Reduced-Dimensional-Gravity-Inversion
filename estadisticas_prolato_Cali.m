%*******************************************************************************
% CLEAR LAKE DATASET - STATISTICAL ANALYSIS
% PROLATE ELLIPSOIDS
%
% Purpose:
%   Performs the statistical analysis of repeated gravity-inversion
%   experiments on the Clear Lake dataset using prolate ellipsoids.
%
% The script computes:
%   - Statistical analysis of inversion error for different model sizes
%   - Boxplots and error histograms
%   - Residual map of the selected best model
%   - Geometrical statistics of the recovered bodies
%   - Median inversion error and IQR versus number of ellipsoids
%
% Authors:
%   GPINV Research Group
%   University of Oviedo, Spain
%
% Copyright (c) 2026 GPINV Research Group, University of Oviedo
%*******************************************************************************

clear
clc
close all

%% ================================================================
% 1. READ ALL RESULTS FILES
% ================================================================

files = dir('Cali_Elip_*.mat');

errores = [];
ncuerpos_all = [];
nombres = cell(length(files),1);

for k = 1:length(files)

    load(files(k).name)

    % Número de cuerpos
    ncuerpos = size(cuerpo_invertido,1);

    % Error final
    err = 100 * min(results.error_hist);

    % Guardar
    errores(end+1) = err;
    ncuerpos_all(end+1) = ncuerpos;
    nombres{k} = files(k).name;

end

grupos = unique(ncuerpos_all);

%% ================================================================
% 2. INVERSION ERROR STATISTICS
% ================================================================

disp(' ')
disp('======================================================')
disp('INVERSION ERROR STATISTICS')
disp('======================================================')

for i = 1:length(grupos)

    n = grupos(i);
    idx = ncuerpos_all == n;
    e = errores(idx);

    fprintf('\n=== %d cuerpos ===\n', n)
    fprintf('Media:   %.3f %%\n', mean(e))
    fprintf('Mediana: %.3f %%\n', median(e))
    fprintf('Std:     %.3f %%\n', std(e))
    fprintf('Min:     %.3f %%\n', min(e))
    fprintf('Max:     %.3f %%\n', max(e))
    fprintf('IQR:     %.3f %%\n', iqr(e))

end

%% ================================================================
% 3. BOXPLOT
% ================================================================

figure

boxplot(errores, ncuerpos_all)

xlabel('Number of bodies')
ylabel('Error (%)')
title('Clear Lake dataset - Prolate Ellipsoids')

grid on

%% ================================================================
% 4. HISTOGRAMS
% ================================================================

figure('Name', 'Histogramas por cuerpo', ...
       'Position', [100 100 1200 600])

n_grupos = length(grupos);
cols = ceil(sqrt(n_grupos));
rows = ceil(n_grupos / cols);

for i = 1:n_grupos

    subplot(rows, cols, i)

    n = grupos(i);
    idx = ncuerpos_all == n;

    histogram(errores(idx), 15, 'FaceColor', '#0072BD')

    title([num2str(n) ' cuerpos'])
    xlabel('Error (%)')
    ylabel('Frecuencia')

    grid on

end

%% ================================================================
% 5. BEST MODEL SELECTION
% ================================================================

target = 16;

idx_target = find(ncuerpos_all == target);

if isempty(idx_target)

    disp(['No se encontraron modelos con ', ...
          num2str(target), ' cuerpos.']);

    return

end

errores_target = errores(idx_target);

[~, idx_best_local] = min(errores_target);

idx_best_global = idx_target(idx_best_local);

best_file = nombres{idx_best_global};

disp(' ')
disp('======================================================')
disp('SELECTED BEST PROLATE-ELLIPSOID MODEL')
disp('======================================================')

[~, name, ext] = fileparts(best_file);
disp([name, ext])

%% ================================================================
% 6. CARGA DEL MEJOR MODELO
% ================================================================

load(best_file);

%% ================================================================
% 7. RESIDUAL MAP
% ================================================================

residual = data.gobs - results.gpre;

figure

scatter(opfun.x, opfun.y, 60, residual, 'filled')
colorbar

xlabel('X')
ylabel('Y')
title('Residual Map (\mu Gal) - Prolate Ellipsoids', 'Interpreter','none')

axis equal
grid on

% figure('Name', 'Mapa Residual')
% 
% scatter(opfun.x, opfun.y, 60, residual, 'filled')
% 
% colormap(jet)
% colorbar
% 
% xlabel('X (m)')
% ylabel('Y (m)')
% 
% title(['Residual Map (\mu Gal) - Best Prolate Model (' , num2str(target) , ' bodies)'], ...
%     'Interpreter', 'none')
% 
% axis equal tight
% grid on

%% ================================================================
% 8. GEOMETRICAL STATISTICS
% ================================================================

disp(' ')
disp('======================================================')
disp('GEOMETRICAL STATISTICS')
disp('======================================================')

z = abs(cuerpo_invertido(:,3));
a = cuerpo_invertido(:,4);
b = cuerpo_invertido(:,5);
theta = cuerpo_invertido(:,6);
phi = cuerpo_invertido(:,7);
rho = cuerpo_invertido(:,8);

% Volumen del elipsoide prolato
V = (4/3) * pi .* a .* (b.^2);

% Masa anómala
M = rho .* V;

fprintf('\nMean depth: %.2f m\n', mean(z))
fprintf('Maximum depth: %.2f m\n', max(z))
fprintf('Minimum depth: %.2f m\n', min(z))
fprintf('Depth standard deviation: %.2f m\n', std(z))

fprintf('\nMean semi-major axis (a): %.2f m\n', mean(a))
fprintf('Mean semi-minor axis (b = c): %.2f m\n', mean(b))

fprintf('\nTotal volume: %.2e m^3\n', sum(V))
fprintf('Mean volume: %.2e m^3\n', mean(V))

fprintf('\nMean density contrast: %.2f kg/m^3\n', mean(rho))

fprintf('\nMean azimuth (theta): %.2f°\n', mean(theta*180/pi))
fprintf('Mean dip angle (phi): %.2f°\n', mean(phi*180/pi))

fprintf('\nTotal anomalous mass: %.2e kg\n', sum(M))

disp(' ')
disp('======================================================')
disp('ANALYSIS COMPLETED')
disp('======================================================')

%% ================================================================
% 9. MEDIAN INVERSION ERROR AND IQR VS NUMBER OF BODIES
% ================================================================

medianas = zeros(size(grupos));
iqrs     = zeros(size(grupos));
means    = zeros(size(grupos));

for i = 1:length(grupos)

    n = grupos(i);
    idx = ncuerpos_all == n;
    e = errores(idx);

    means(i)    = mean(e);
    medianas(i) = median(e);
    iqrs(i)     = iqr(e);

end

figure
hold on
box on

yyaxis left

plot(grupos, medianas, '-o', 'LineWidth', 2.2, 'MarkerSize', 7)

ylabel('Median  (%)')

yyaxis right

plot(grupos, iqrs, '-s', 'LineWidth', 2.2, 'MarkerSize', 7)

ylabel('IQR (%)')

xlabel('Number of prolate ellipsoids')

grid on
set(gca,'FontSize',12)

legend('Median error', 'IQR', 'Location', 'northeast')