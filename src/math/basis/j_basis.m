function J = j_basis(L, x)
%J_BASIS Evaluate the Fourier-Legendre basis functions.
%
%   J = j_basis(L, x)
%
%   Inputs
%   ------
%   L : scalar
%       Number of basis functions.
%
%   x : array
%       Evaluation points.
%
%   Output
%   ------
%   J : L x numel(x) array
%       Values of the Fourier-Legendre basis functions evaluated at x.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        L
        x
    end


    %% Initialize basis matrix

    J = zeros(L, numel(x));
    ax = abs(x);


    %% Evaluate basis functions

    for n = 0:L-1
        base = real(besselj(n + 0.5, 2*pi*ax) .* sqrt(1 ./ (4*ax)));

        base(x == 0) = double(n == 0);
        base(x < 0) = (-1)^n * base(x < 0);

        J(n + 1, :) = (-1i)^n * base;
    end

end