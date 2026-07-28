%*******************************************************************************
% SYNTHETIC EXPERIMENTS - STATISTICAL ANALYSIS AND RESIDUAL MAPS
%
% Purpose:
%   Performs the statistical analysis of the synthetic gravity-inversion
%   experiments for the three geometric parametrizations considered in the
%   associated manuscript.
%
% The script computes:
%   - Mean and median inversion error
%   - Standard deviation and interquartile range (IQR)
%   - Best inversion error
%   - Average CPU time
%   - Boxplots and error histograms
%   - Residual maps for the best model of each parametrization
%
% Residual definition:
%   residual = gobs - gpre
%
% Parametrizations:
%   1 - Prolate ellipsoids
%   2 - Square-based prisms
%   3 - Rectangular prisms
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
% FILE PREFIXES
% ================================================================

prefix_elip = 'SYN_Elip_*.mat';
prefix_sq   = 'SYN_SqPrism_*.mat';
prefix_rect = 'SYN_Pris_rect_*.mat';

labels = {'Prolate Ellipsoids', ...
          'Square Prisms', ...
          'Rectangular Prisms'};

prefixes = {prefix_elip, prefix_sq, prefix_rect};

%% ================================================================
% STORAGE
% ================================================================

all_errors = cell(1,3);
all_times  = cell(1,3);

fprintf('\n====================================================\n')
fprintf('SYNTHETIC EXPERIMENTS - STATISTICAL ANALYSIS\n')
fprintf('====================================================\n')

%% ================================================================
% LOOP OVER GEOMETRIES
% ================================================================

for g = 1:3
    
    files = dir(prefixes{g});
    
    errores = [];
    tiempos = [];
    
    fprintf('\n--------------------------------------------\n')
    fprintf('%s\n', labels{g})
    fprintf('--------------------------------------------\n')
    
    for k = 1:length(files)
        
        load(files(k).name)
        
        %----------------------------------------------------------
        % % Final inversion error (%)
        %----------------------------------------------------------
        err = 100 * min(results.error_hist);
        errores(end+1) = err;
        
        %----------------------------------------------------------
        % CPU time
        %----------------------------------------------------------
        if exist('Duracion','var')
            tiempos(end+1) = Duracion;
        elseif isfield(results,'time')
            tiempos(end+1) = results.time;
        else
            tiempos(end+1) = NaN;
        end
        
        clear results Duracion data gpre cuerpo_invertido
    end
    
    all_errors{g} = errores;
    all_times{g}  = tiempos;
    
    fprintf('Number of runs: %d\n', length(errores))
    fprintf('Mean error: %.4f\n', mean(errores))
    fprintf('Median error: %.4f\n', median(errores))
    fprintf('Std: %.4f\n', std(errores))
    fprintf('Min (Best error): %.4f\n', min(errores))
    fprintf('Max: %.4f\n', max(errores))
    fprintf('IQR: %.4f\n', iqr(errores))
    fprintf('Average CPU Time: %.4f s\n', mean(tiempos,'omitnan'))
    
end

%% ================================================================
% BOXPLOT COMPARISON
% ================================================================

figure

errors_concat = [all_errors{1}, all_errors{2}, all_errors{3}];

group = [ones(1,length(all_errors{1})), ...
         2*ones(1,length(all_errors{2})), ...
         3*ones(1,length(all_errors{3}))];

boxplot(errors_concat, group, 'Labels', labels)

ylabel('Error (%)')
title('Synthetic Experiments')
grid on

%% ================================================================
% HISTOGRAMS
% ================================================================

figure

for g = 1:3
    
    subplot(1,3,g)
    
    histogram(all_errors{g}, 15)
    
    title(labels{g})
    xlabel('Error (%)')
    ylabel('Frequency')
    grid on
    
end

sgtitle('Synthetic Experiments')

%% ================================================================
% RESIDUAL MAPS - BEST MODEL
% ================================================================
%
% For each parametrization:
%   - Find best run (minimum inversion error)
%   - Recompute forward response from best model
%   - Compute residual:
%
%       residual = gobs - gpre_best
%
% Using:
%   best model = cuerpo_invertido
%
% This avoids inconsistencies caused by stored gpre values.
%
%% ================================================================

figure('Name','Residual Map - Best Model', ...
       'Position',[100 100 1500 500])

sgtitle('Residual Maps - Best Models')

%% ------------------------------------------------
% FIRST PASS:
% Determine common symmetric color scale
%% ------------------------------------------------

global_residuals = [];

for g = 1:3
    
    files = dir(prefixes{g});
    
    best_error = inf;
    best_residual = [];
    
    for k = 1:length(files)
        
        load(files(k).name)
        
        err = 100 * min(results.error_hist);
        
        if err < best_error
            
            best_error = err;
            
            if exist('data','var') && ...
               isfield(data,'gobs') && ...
               exist('cuerpo_invertido','var')
           
                %----------------------------------
                % Recompute predicted anomaly
                %----------------------------------
                
                if g == 1
                    % Prolate Ellipsoids
                    
                    [gpre_temp, ~] = ...
                        GravElipProlato(cuerpo_invertido, ptos);
                    
                    gpre_best = -gpre_temp(:,3) * 1.0e8;
                    
                elseif g == 2
                    % Square Prisms
                    
                    [gpre_temp, ~, ~, ~] = ...
                        GravPrismaCuadrado(cuerpo_invertido, ...
                                           ptos, 1, 1.0);
                    
                    gpre_best = -gpre_temp(:,3) * 1.0e8;
                    
                elseif g == 3
                    % Rectangular Prisms
                    
                    [gpre_temp, ~, ~, ~] = ...
                        GravPrismaRectangular(cuerpo_invertido, ...
                                              ptos, 1, 1.0);
                    
                    gpre_best = -gpre_temp(:,3) * 1.0e8;
                    
                end
                
                %----------------------------------
                % True residual
                %----------------------------------
                
                best_residual = data.gobs - gpre_best;
                best_residual = best_residual(:);
                
            end
        end
        
        clear results Duracion data gpre gpre_best cuerpo_invertido
    end
    
    if ~isempty(best_residual)
        global_residuals = [global_residuals; best_residual];
    end
    
end

% Symmetric color scale

cmax = max(abs(global_residuals));

%% ------------------------------------------------
% SECOND PASS:
% Plot best residual map for each geometry
%% ------------------------------------------------

for g = 1:3
    
    files = dir(prefixes{g});
    
    best_error = inf;
    best_residual = [];
    best_X = [];
    best_Y = [];
    
    for k = 1:length(files)
        
        load(files(k).name)
        
        err = 100 * min(results.error_hist);
        
        if err < best_error
            
            best_error = err;
            
            if exist('data','var') && ...
               isfield(data,'gobs') && ...
               exist('cuerpo_invertido','var')
           
                %----------------------------------
                % Recompute predicted anomaly
                %----------------------------------
                
                if g == 1
                    % Prolate Ellipsoids
                    
                    [gpre_temp, ~] = ...
                        GravElipProlato(cuerpo_invertido, ptos);
                    
                    gpre_best = -gpre_temp(:,3) * 1.0e8;
                    
                elseif g == 2
                    % Square Prisms
                    
                    [gpre_temp, ~, ~, ~] = ...
                        GravPrismaCuadrado(cuerpo_invertido, ...
                                           ptos, 1, 1.0);
                    
                    gpre_best = -gpre_temp(:,3) * 1.0e8;
                    
                elseif g == 3
                    % Rectangular Prisms
                    
                    [gpre_temp, ~, ~, ~] = ...
                        GravPrismaRectangular(cuerpo_invertido, ...
                                              ptos, 1, 1.0);
                    
                    gpre_best = -gpre_temp(:,3) * 1.0e8;
                    
                end
                
                %----------------------------------
                % Residual
                %----------------------------------
                
                best_residual = data.gobs - gpre_best;
                best_residual = best_residual(:);
                
                % Coordinates
                
                best_X = ptos(:,1);
                best_Y = ptos(:,2);
                
            end
        end
        
        clear results Duracion data gpre gpre_best cuerpo_invertido
    end
    
    %% --------------------------------------------
    % Plot
    %% --------------------------------------------
    
    subplot(1,3,g)
    
    scatter(best_X, best_Y, ...
            90, best_residual, ...
            'filled', ...
            'MarkerEdgeColor','k')
    
    title(sprintf('%s\nBest inversion error = %.2f%%', ...
      labels{g}, best_error), ...
      'FontWeight','bold')
    
    xlabel('X')
    ylabel('Y')
    
    axis equal
    grid on
    box on
    
    colormap(jet)
    caxis([-cmax cmax])
    
    colorbar
    
end

%% Optional
% saveas(gcf,'residual_map_best_model.png')

%% ================================================================
% LATEX TABLE VALUES
% ================================================================

fprintf('\n====================================================\n')
fprintf('VALUES FOR LATEX TABLE\n')
fprintf('====================================================\n')

fprintf('\nCopy these values into your paper table:\n\n')

for g = 1:3
    
    e = all_errors{g};
    t = all_times{g};
    
    fprintf('%s\n', labels{g})
    fprintf('Median error = %.3f\n', median(e))
    fprintf('IQR         = %.3f\n', iqr(e))
    fprintf('Best error   = %.3f\n', min(e))
    fprintf('CPU Time    = %.3f s\n', mean(t,'omitnan'))
    fprintf('\n')
    
end