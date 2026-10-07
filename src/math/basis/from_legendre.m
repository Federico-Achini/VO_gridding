function fun = from_legendre(c, N, endpoints)
%FROM_LEGENDRE Reconstruct a function from Legendre coefficients.
%
%   fun = from_legendre(c, N, endpoints)
%
%   Inputs
%   ------
%   c : array
%       Legendre coefficients. The first column contains the real
%       coefficients and the second column contains the imaginary
%       coefficients.
%
%   N : scalar
%       Number of points at which the function is evaluated.
%
%   endpoints : character vector
%       Endpoint convention used for the discretization. Default: 'left'.
%
%   Output
%   ------
%   fun : array
%       Reconstructed function values.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        c
        N
        endpoints = 'left'
    end


    %% Legendre basis

    [L, ~] = size(c);
    tau = discretize(-1, 1, N, endpoints);
    M = legendre_poly(L, tau).';


    %% Function reconstruction

    fun = M * (c(:, 1) + 1i*c(:, 2));

end