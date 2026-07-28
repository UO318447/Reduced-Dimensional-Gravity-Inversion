%% QuickTest.m
%==========================================================================
% Quick test for the Reduced-Dimensional Gravity Inversion software
%
% This script performs a short RR-GPSO inversion to verify that the
% repository has been correctly installed.
%
% Expected running time:
%   approximately 1–2 minutes on a standard desktop computer.
%==========================================================================

clear
clc
close all

addpath(genpath(pwd))

fprintf('\n');
fprintf('=============================================\n');
fprintf(' Reduced-Dimensional Gravity Inversion\n');
fprintf(' Quick Test\n');
fprintf('=============================================\n\n');

% Short execution parameters
n_runs        = 1;
noise_level   = 0.10;
quick_maxiter = 20;
quick_swarm   = 20;

MAIN_PSO_SINTETICO

fprintf('\n');
fprintf('Quick test finished successfully.\n');
fprintf('The software is correctly installed.\n');