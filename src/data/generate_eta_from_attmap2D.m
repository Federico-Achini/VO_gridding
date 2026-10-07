function eta = generate_eta_from_attmap2D(a_map, x, y, k, chi)
%GENERATE_ETA_FROM_ATTMAP2D Generate a target error profile from a 2-D
%attention map.
%
% Author: Federico Achini
% Date:   August 2025
%
%   eta = generate_eta_from_attmap2D(a_map, x, y, k, chi)
%
%   Inputs
%   ------
%   a_map : 2-D array
%       Attention map.
%
%   x : vector
%       Coordinates along the first dimension.
%
%   y : vector
%       Coordinates along the second dimension.
%
%   k : scalar
%       Scaling parameter controlling the target error.
%
%   chi : scalar
%       Exponent applied to the normalized attention map.
%
%   Output
%   ------
%   eta : 2-D array
%       Target error profile defined on the x-y grid.
%
%   ----------------------------------------------------------------------

    %% Grid dimensions

    Mx = length(x);
    My = length(y);


    %% Resize attention map

    a_map = imresize(a_map, [Mx, My]);


    %% Normalize attention map

    norms = sum(abs(a_map), "all");

    a_map = ...
        (a_map ./ norms * Mx * My).^chi;


    %% Generate target error profile

    eta = 10.^(-k * a_map);

end