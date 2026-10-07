function P = legendre_poly(L, x)
%LEGENDRE_POLY Evaluate the first L Legendre polynomials P_0..P_{L-1}.
%
%   P = legendre_poly(L, x)
%
%   Inputs
%   ------
%   L : positive integer
%       Number of Legendre polynomials to evaluate.
%
%   x : vector
%       Evaluation points.
%
%   Output
%   ------
%   P : L x numel(x) array
%       Legendre polynomial values, with
%       P(n+1, :) = P_n(x).
%
%   Description
%   -----------
%   The Legendre polynomials are evaluated using the stable three-term
%   recurrence
%
%       (n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x).
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        L (1,1) double {mustBePositive, mustBeInteger}
        x double
    end


    %% Initialize evaluation grid

    x = x(:).';              % 1 x numel(x) row vector
    P = zeros(L, numel(x));

    P(1, :) = 1;             % P_0(x) = 1

    if L > 1
        P(2, :) = x;         % P_1(x) = x
    end


    %% Three-term recurrence

    for n = 1:L-2
        P(n+2, :) = ((2*n+1)*x.*P(n+1, :) - n*P(n, :)) / (n+1);
    end

end