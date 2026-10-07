function eta = generate_eta_from_attmap(a_map, x, nScan, k, chi)
%GENERATE_ETA_FROM_ATTMAP Generate a target error profile from an attention map.
%
% Author: Federico Achini
% Date:   August 2025
%
%   eta = generate_eta_from_attmap(a_map, x, nScan, k, chi)
%
%   Inputs
%   ------
%   a_map : 2-D array
%       Attention map.
%
%   x : vector
%       Coordinates along the first dimension.
%
%   nScan : scalar
%       Number of scan positions along the second dimension.
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
%       Target error profile.
%
%   ----------------------------------------------------------------------

    %% Grid dimension

    M = length(x);


    %% Resize attention map

    a_map = imresize(a_map, [M, nScan]);


    %% Normalize attention map along each scan

    norms = sum(abs(a_map), 1);

    a_map = (a_map ./ norms).^chi * M;


    %% Generate target error profile

    eta = 10.^(-k .* a_map);

end