function coord = radial_trajectory(ntr, npts)
%RADIAL_TRAJECTORY Generate a 2-D radial k-space sampling trajectory.
%
%   coord = radial_trajectory(ntr, npts)
%
%   Inputs
%   ------
%   ntr : scalar
%       Number of radial spokes (trajectories).
%
%   npts : scalar
%       Number of sample points per spoke.
%
%   Output
%   ------
%   coord : 2 x npts x ntr array
%       Coordinates of the radial trajectory:
%           coord(1, :, :) contains the y-components.
%           coord(2, :, :) contains the x-components.
%
%   Description
%   -----------
%   Generates a polar radial trajectory in 2-D k-space. The angles are
%   distributed according to the golden-ratio increment, and the radial
%   coordinate is generated over the specified sampling interval.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Input arguments

    arguments
        ntr   % Number of radial spokes
        npts  % Number of points per spoke
    end


    %% Trajectory parameters

    phi = pi * (1 + sqrt(5));           % Using golden ratio
    theta = phi * (0:ntr-1);            % Uniformly spaced angular increment
    r = discretize(0, 0.5, npts);       % Radial sampling points


    %% Allocate trajectory

    coord = zeros(2, npts, ntr);


    %% Compute Cartesian coordinates

    % Compute x and y coordinates using polar-to-Cartesian conversion.
    coord(1, :, :) = reshape(r .* cos(theta), 1, npts, ntr);   % x-coordinates
    coord(2, :, :) = reshape(r .* sin(theta), 1, npts, ntr);   % y-coordinates

end