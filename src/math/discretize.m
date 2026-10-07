function y = discretize(a, b, N, endpoint)
%DISCRETIZE Uniformly discretize an interval into N points.
%
%   y = discretize(a, b, N)
%
%   Inputs
%   ------
%   a : scalar
%       Start of the interval.
%
%   b : scalar
%       End of the interval.
%
%   N : scalar
%       Number of points to discretize the interval into.
%
%   endpoint : string or character vector
%       Endpoint convention:
%           'left'  - include a, exclude b
%           'right' - exclude a, include b
%           'none'  - exclude both endpoints
%           'both'  - include both endpoints
%
%       Default: "left".
%
%   Output
%   ------
%   y : column vector
%       Uniformly spaced discretization points.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        a
        b
        N
        endpoint = "left"
    end


    %% Discretization

    % Remove the last point to exclude endpoint b (for periodicity).
    if endpoint == "left"

        y = linspace(a, b, N + 1);
        y = y(1:end-1);

    elseif endpoint == "right"

        y = linspace(a, b, N + 1);
        y = y(2:end);

    elseif endpoint == "none"

        y = linspace(a, b, N + 2);
        y = y(2:end-1);

    elseif endpoint == "both"

        y = linspace(a, b, N);

    else

        error("Choose an endpoint option.")

    end


    %% Output orientation

    % Convert row vector to column vector.
    y = y.';

end