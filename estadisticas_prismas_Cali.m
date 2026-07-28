%*******************************************************************************
% CLEAR LAKE DATASET - STATISTICAL ANALYSIS
% RECTANGULAR PRISMS
%
% Purpose:
%   Performs the statistical analysis of repeated gravity-inversion
%   experiments on the Clear Lake dataset using rectangular prisms.
%
% The script computes:
%   - Statistical analysis of inversion error for different model sizes
%   - Boxplots and error histograms
%   - Residual map of the selected best model
%   - Geometrical statistics of the recovered bodies
%   - Median inversion error and IQR versus number of prisms
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

files = dir('Cali_PRect_*.mat');

errores = [];
ncuerpos_all = [];
nombres = {};

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
disp('ESTADÍSTICAS DEL ERROR')
disp('======================================================')

for i = 1:length(grupos)

    n = grupos(i);
    idx = ncuerpos_all == n;
    e = errores(idx);

    fprintf('\n=== %d cuerpos ===\n', n)
    fprintf('Media: %.3f\n', mean(e))
    fprintf('Mediana: %.3f\n', median(e))
    fprintf('Std: %.3f\n', std(e))
    fprintf('Min: %.3f\n', min(e))
    fprintf('Max: %.3f\n', max(e))
    fprintf('IQR: %.3f\n', iqr(e))

end

%% ================================================================
% 3. BOXPLOT
% ================================================================

figure

boxplot(errores, ncuerpos_all)

xlabel('Number of bodies')
ylabel('Error (%)')
%title('Distribución estadística del error')
title('Clear Lake dataset - Rectangular Prisms')
grid on

%% ================================================================
% 4. HISTOGRAMS
% ================================================================

figure

for i = 1:length(grupos)

    subplot(2,2,i)

    n = grupos(i);
    idx = ncuerpos_all == n;

    histogram(errores(idx),15)

    title([num2str(n) ' cuerpos'])
    xlabel('Error (%)')
    ylabel('Frequency')
    grid on

end

%% ================================================================
% 5. BEST MODEL SELECTION
% ================================================================

target = 4;

idx_target = find(ncuerpos_all == target);

errores_target = errores(idx_target);

[~, idx_best_local] = min(errores_target);

idx_best_global = idx_target(idx_best_local);

best_file = nombres{idx_best_global};

disp(' ')
disp('======================================================')
disp('BEST MODEL SELECTION')
disp('======================================================')
disp(best_file)

%% ================================================================
% 6. LOAD BEST MODEL
% ================================================================

load(best_file)

%% ================================================================
% 7. RESIDUAL MAP
% ================================================================

residual = data.gobs - results.gpre;

figure

scatter(opfun.x, opfun.y, 60, residual, 'filled')
colorbar

xlabel('X')
ylabel('Y')
title('Residual Map - Best Model', 'Interpreter','none')

axis equal
grid on

%% ================================================================
% 8. GEOMETRICAL STATISTICS
% ================================================================

disp(' ')
disp('======================================================')
disp('GEOMETRICAL STATISTICS')
disp('======================================================')

% Parámetros del prisma:
% [x y z A B C theta phi rho]

z = abs(cuerpo_invertido(:,3));

A = cuerpo_invertido(:,4);
B = cuerpo_invertido(:,5);
C = cuerpo_invertido(:,6);

theta = cuerpo_invertido(:,7);
phi = cuerpo_invertido(:,8);

rho = cuerpo_invertido(:,9);

%% Volumen

V = A .* B .* C;

%% Masa anómala

M = rho .* V;

%% Resultados

fprintf('\nMean depth: %.2f m\n', mean(z))
fprintf('Maximum depth: %.2f m\n', max(z))
fprintf('Minimum depth: %.2f m\n', min(z))
fprintf('Depth standard deviation: %.2f m\n', std(z))

fprintf('\nTotal volume: %.2e m^3\n', sum(V))
fprintf('Mean volume: %.2e m^3\n', mean(V))

fprintf('\nMean density contrast: %.2f kg/m^3\n', mean(rho))

fprintf('\nMean azimuth (theta): %.2f°\n', mean(theta*180/pi))
fprintf('Mean dip angle (phi): %.2f°\n', mean(phi*180/pi))

fprintf('\nTotal anomalous mass: %.2e kg\n', sum(M))

disp(' ')
disp('======================================================')
disp('ANÁLISIS COMPLETADO')
disp('======================================================')


%% ================================================================
%  FIGURA RMSE + IQR vs Nº DE CUERPOS
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
plot(grupos, medianas, '-o', ...
    'LineWidth', 2.2, ...
    'MarkerSize', 7)

ylabel('Median error (%)')

yyaxis right
plot(grupos, iqrs, '-s', ...
    'LineWidth', 2.2, ...
    'MarkerSize', 7)

ylabel('IQR (%)')

xlabel('Number of rectangular prisms')

%title('Sensitivity analysis: inversion quality vs model complexity')

grid on
set(gca,'FontSize',12)

legend('Median error','IQR','Location','northeast')