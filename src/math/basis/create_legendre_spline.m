function spl = create_legendre_spline(c, W, D)
%CREATE_LEGENDRE_SPLINE Create a spline representation from Legendre
%coefficients.
%
%   spl = create_legendre_spline(c, W, D)
%
%   Inputs
%   ------
%   c : array
%       Legendre coefficients.
%
%   W : scalar
%       Kernel half-width.
%
%   D : scalar
%       Discretization factor.
%
%   Output
%   ------
%   spl : structure or object
%       Spline representation of the kernel.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Discretization

    N = 2*W*D + 1;
    breaks = discretize(-W, W, N, 'both');


    %% Legendre reconstruction

    y = from_legendre(c, N, 'both');


    %% Spline representation

    spl = real_imag_spline(breaks, y);

end