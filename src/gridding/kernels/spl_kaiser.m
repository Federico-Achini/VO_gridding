function spl = spl_kaiser(W, D)
%SPL_KAISER Construct a spline representation of a Kaiser-Bessel kernel.
%
%   spl = spl_kaiser(W, D)
%
%   Inputs
%   ------
%   W : scalar
%       Kernel width.
%
%   D : scalar
%       Number of samples per unit interval.
%
%   Output
%   ------
%   spl : structure
%       Spline representation of the Kaiser-Bessel kernel.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Kaiser-Bessel parameters

    m = 0;
    beta = 4.68*W;


    %% Kernel grid

    x = discretize(-W, W, 2*W*D + 1, 'both');
    f = sqrt(1 - (x/W).^2);


    %% Kaiser-Bessel kernel

    denom = besseli(m, beta);
    kb = f.^m .* besseli(m, beta*f) / denom;


    %% Spline representation

    spl = real_imag_spline(x, real(kb));

end