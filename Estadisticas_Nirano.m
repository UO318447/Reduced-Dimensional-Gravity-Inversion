%*******************************************************************************
% NIRANO DATASET - STATISTICAL ANALYSIS
%
% Purpose:
%   Performs the statistical analysis of repeated gravity-inversion
%   experiments on the Nirano dataset using two geometric parametrizations.
%
% The script computes:
%   - Mean and median inversion error
%   - Standard deviation and interquartile range (IQR)
%   - Best inversion error
%   - Boxplots and histograms
%   - Residual maps of the best models
%   - Wilcoxon rank-sum statistical test
%
% Parametrizations:
%   1 - Rectangular prisms
%   2 - Prolate ellipsoids
%
% Authors:
%   GPINV Research Group
%   University of Oviedo, Spain
%
% Copyright (c) 2026 GPINV Research Group, University of Oviedo
%*******************************************************************************

%% ================================================================
% PART A - RECTANGULAR PRISMS
% ================================================================

clear
clc
close all

files = dir('rn_Nirano_PRect_*.mat');

errores = [];
nombres = {};

for k = 1:length(files)

    load(files(k).name)

    err = 100 * min(results.error_hist);

    errores(end+1) = err;
    nombres{k} = files(k).name;

end

fprintf('\n===============================================\n')
fprintf('NIRANO - RECTANGULAR PRISMS\n')
fprintf('===============================================\n')

fprintf('Número de runs: %d\n', length(errores))
fprintf('Mean: %.3f\n', mean(errores))
fprintf('Median: %.3f\n', median(errores))
fprintf('Std: %.3f\n', std(errores))
fprintf('Min: %.3f\n', min(errores))
fprintf('Max: %.3f\n', max(errores))
fprintf('IQR: %.3f\n', iqr(errores))

errores_prismas = errores;

% figure
% boxplot(errores)
% ylabel('Error (%)')
% title('Nirano - Prismas - Boxplot')
% grid on
% 
% figure
% histogram(errores,15)
% xlabel('Error (%)')
% ylabel('Frecuencia')
% title('Nirano - Prismas - Histograma')
% grid on

[~, idx_best] = min(errores);
best_file = nombres{idx_best};

fprintf('\nMejor modelo prismas: %s\n', best_file)

load(best_file)

residual = data.gobs - results.gpre;

figure
scatter(opfun.x, opfun.y, 70, residual, 'filled')
colorbar
xlabel('X')
ylabel('Y')
title('Residual Map - Best Rectangular-Prism Model')
axis equal
grid on


%% ================================================================
% PART B - PROLATE ELLIPSOIDS
% ================================================================

files = dir('rn_Nirano_Elip_*.mat');

errores = [];
nombres = {};

for k = 1:length(files)

    load(files(k).name)

    err = 100 * min(results.error_hist);

    errores(end+1) = err;
    nombres{k} = files(k).name;

end

fprintf('\n===============================================\n')
fprintf('NIRANO - PROLATE ELLIPSOIDS\n')
fprintf('===============================================\n')

fprintf('Mean inversion error: %.3f %%\n', mean(errores))
fprintf('Median inversion error: %.3f %%\n', median(errores))
fprintf('Standard deviation: %.3f %%\n', std(errores))
fprintf('Best inversion error: %.3f %%\n', min(errores))
fprintf('Worst inversion error: %.3f %%\n', max(errores))
fprintf('Interquartile range (IQR): %.3f %%\n', iqr(errores))
errores_prolatos = errores;

% figure
% boxplot(errores)
% ylabel('Error (%)')
% title('Nirano - Prolatos - Boxplot')
% grid on
% 
% figure
% histogram(errores,15)
% xlabel('Error (%)')
% ylabel('Frecuencia')
% title('Nirano - Prolatos - Histograma')
% grid on

[~, idx_best] = min(errores);
best_file = nombres{idx_best};

fprintf('\nBest rectangular-prism model: %s\n', best_file)

load(best_file)

residual = data.gobs - results.gpre;

figure
scatter(opfun.x, opfun.y, 70, residual, 'filled')
colorbar
xlabel('X')
ylabel('Y')
title('Best model - Prolate Ellipsoids')
axis equal
grid on


%% ================================================================
% BOXPLOT COMPARATIVO FINAL
% ================================================================

figure

errores_total = [errores_prismas(:); errores_prolatos(:)];

grupo = [ ...
    repmat({'Rectangular Prisms'}, length(errores_prismas), 1);...
    repmat({'Prolate Ellipsoids'}, length(errores_prolatos), 1);    
    ];

boxplot(errores_total, grupo)

ylabel('Error (%)')
title('Nirano dataset')
grid on

%% ================================================================
% COMPARACIÓN ESTADÍSTICA ENTRE GEOMETRÍAS
% ================================================================

[p,h,stats] = ranksum(errores_prismas, errores_prolatos);

disp(' ')
disp('======================================================')
disp('WILCOXON RANK-SUM TEST')
disp('======================================================')

fprintf('p-value = %.6f\n', p);

if h
    disp('Result: Statistically significant difference (α = 0.05)')
else
    disp('Result: No statistically significant difference (α = 0.05)')
end

fprintf('Rank-sum statistic = %.2f\n', stats.ranksum);


%% ================================================================
% VALUES FOR THE MANUSCRIPT TABLE
% ================================================================
%
% These values can be directly incorporated into the statistical
% comparison table reported in the manuscript
%
% -------------------------------------------------------------
% Geometry | Mean | Median | Std | IQR | Best Error
% -------------------------------------------------------------
% Prismas  |  ... |   ...  | ... | ... |    ...
% Prolatos |  ... |   ...  | ... | ... |    ...
% -------------------------------------------------------------
%

