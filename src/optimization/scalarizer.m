function y = scalarizer(x, p, rho, eps_reg)
%SCALARIZER Compute a scalar objective value from a vector based on norms.
%
%   y = scalarizer(x, p, rho)
%   y = scalarizer(x, p, rho, eps_reg)
%
%   Inputs
%   ------
%   x : vector
%       Input vector of values (e.g., errors).
%
%   p : scalar
%       Norm order (finite positive scalar, e.g., 1, 2).
%
%   rho : scalar
%       Scalarization factor controlling positive-part weighting.
%
%   eps_reg : scalar, optional
%       Regularization floor used only when p <= 1.
%       Default: 1e-25.
%
%   Output
%   ------
%   y : scalar
%       Scalar objective value combining full and positive parts of x.
%
%   Description
%   -----------
%   This function calculates a scalar value by combining the L^p norm
%   of the full vector x and the L^p norm of the positive elements of x.
%   The positive-part norm is weighted by rho and subtracted from the
%   full norm.
%
%   For p <= 1, |x_i|^p (and x_pos_i^p) has an unbounded derivative as
%   x_i -> 0, which causes gradient/Hessian blow-up near-zero residuals
%   -- exactly the residuals this scalarizer is often applied to. To
%   keep the value finite and its derivative bounded there, each term is
%   regularized as (|x_i| + eps_reg)^p - eps_reg^p, which matches
%   |x_i|^p closely for |x_i| >> eps_reg and vanishes smoothly at
%   x_i = 0 instead of having a singular slope/curvature.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Regularization

    if nargin < 4 || isempty(eps_reg)

        eps_reg = 1e-25;   % Only used when p <= 1; tune relative to your
                           % expected residual floor (see note below)

    end


    %% Positive part

    x_pos = x(x > 0);


    %% Compute p-power terms

    if p <= 1

        % Regularized quasi-norm terms: bounded derivative as elements -> 0.
        normp_full = sum((abs(x) + eps_reg).^p - eps_reg^p);
        normp_pos  = sum((x_pos + eps_reg).^p - eps_reg^p);

    else

        % p > 1: |x|^p is already smooth at 0, no regularization needed.
        normp_full = norm(x, p)^p;
        normp_pos  = norm(x_pos, p)^p;

    end


    %% Scalarized objective

    %   -( ||x||_p^p - rho * ||x_pos||_p^p )
    if rho == inf

        y = normp_pos;

    else

        y = -(normp_full - rho * normp_pos);

    end

end