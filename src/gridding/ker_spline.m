function spl = ker_spline(c, W, D, mode)
%KER_SPLINE Convert kernel coefficients to a spline representation.
%
%   spl = ker_spline(c, W, D)
%
%   Inputs
%   ------
%   c : array
%       Kernel coefficients.
%
%   W : scalar
%       Half-bandwidth parameter defining the support interval.
%
%   D : scalar
%       Discretization factor.
%
%   mode : string
%       Kernel representation mode. Default: "slepian".
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

    arguments
        c
        W
        D
        mode = "slepian"
    end


    %% Discretization

    % Total number of discretization points in the interval [-W, W).
    T = 2 * W * D;

    % Uniform discretization of interval [-W, W).
    nu = discretize(-W, W, T + 1, 'both');


    %% Kernel reconstruction

    if mode == "slepian"

        % Reconstruct kernel function by weighted sum of Slepian sequences.
        ker = slepian(c, T + 1, W);

        % Generate spline representations for real and imaginary parts.
        spl = real_imag_spline(nu, ker);

    elseif mode == "vector"

        c = c(:, 1) + 1i * c(:, 2);
        nu = [nu; W];

        spl = real_imag_spline(nu, c);

    end

end