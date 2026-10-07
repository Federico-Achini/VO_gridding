% --------------------
% |  my_startup.m    |
% --------------------
%
% This script initializes the MATLAB session and loads the global
% plotting and visualization defaults.
%
% It closes open figures, clears the workspace, and clears the
% command window before running the custom startup configuration.
%
% Author: Federico Achini
% Date:   August 2025
%
% ------------------------------------------------------------------------

%% Initialize MATLAB session

% Close all open figures.
close all;

% Clear all variables from the workspace.
clear all;

% Clear the command window.
clc;


%% Configure plotting defaults

% Load global plotting and visualization settings.
run my_startup();