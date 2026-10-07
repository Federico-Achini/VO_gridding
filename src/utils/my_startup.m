% -------------------
% |  my_startup.m   |
% -------------------
%
% This script configures global MATLAB plotting defaults to ensure
% consistent, publication-quality figure appearance.
%
% Settings include LaTeX text rendering, font sizes, grid display,
% line and marker styles, minor ticks/grids, and window style.
%
% It also defines helper functions to retrieve custom blue and red
% color codes by index for use in plotting.
%
% Author: Federico Achini
% Date:   August 2025
%
% ------------------------------------------------------------------------

%% Text rendering

% Use LaTeX interpreter for all text objects (titles, labels, legends, etc.).
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');


%% Default font sizes

set(groot, ...
    'defaultAxesFontSize', 16, ...          % Font size for axes labels and ticks
    'defaultTextFontSize', 20, ...           % Font size for standalone text objects
    'defaultTextFontSizeMode', 'manual' ...  % Fix font size, disable auto-scaling
);


%% Grid settings

% Enable grid lines for X and Y axes by default.
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');

% Enable minor ticks on both X and Y axes.
set(groot, 'DefaultAxesXMinorTick', 'on', ...
          'DefaultAxesYMinorTick', 'on');

% Enable minor grid lines on X axis, with manual control over the mode.
set(groot, 'defaultAxesXMinorGrid', 'on', ...
          'defaultAxesXMinorGridMode', 'manual');


%% Line and marker settings

% Set default line width for plot lines.
set(groot, 'DefaultLineLineWidth', 1.5);

% Set default marker size for plot markers (circle, x, etc.).
set(groot, 'defaultLineMarkerSize', 6);


%% Figure window style

% Set default figure window style to docked (inside MATLAB desktop).
set(0, 'DefaultFigureWindowStyle', 'docked');


%% Color palette handles

% Define function handles for blue, red, and green color palettes.
blue = @(i) m_blue(i);
red   = @(i) m_red(i);
green = @(i) m_green(i);


%% Blue palette

% Function: m_blue
% Purpose: Returns a hex color string for a custom blue shade by index.
%
% Inputs:
%   i - index specifying shade (0: light, 1: medium, 2: dark)
%
% Outputs:
%   y - string with hex color code for the selected blue shade

function y = m_blue(i)

    switch i
        case 0
            y = "#00CCFF";   % Light blue
        case 1
            y = "#0066FF";   % Medium blue
        case 2
            y = "#003399";   % Dark blue
        otherwise
            y = "#0066FF";   % Default medium blue if index invalid
    end

end


%% Red palette

% Function: m_red
% Purpose: Returns a hex color string for a custom red shade by index.
%
% Inputs:
%   i - index specifying shade (0: light/orange, 1: medium, 2: dark/brown)
%
% Outputs:
%   y - string with hex color code for the selected red shade

function y = m_red(i)

    switch i
        case 0
            y = "#FF9966";   % Light red / orange
        case 1
            y = "#CC3300";   % Medium red
        case 2
            y = "#663300";   % Dark red / brownish
        otherwise
            y = "#CC3300";   % Default medium red if index invalid
    end

end


%% Green palette

% Function: m_green
% Purpose: Returns a hex color string for a custom green shade by index.
%
% Inputs:
%   i - index specifying shade (0: light/orange, 1: medium, 2: dark/brown)
%
% Outputs:
%   y - string with hex color code for the selected green shade

function y = m_green(i)

    switch i
        case 0
            y = "#ADFF2F";   % Light green / yellow
        case 1
            y = "#32CD32 ";  % Medium green
        case 2
            y = "#008000 ";  % Dark green
        otherwise
            y = "#32CD32";   % Default medium green if index invalid
    end

end