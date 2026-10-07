function [f, g, H] = precond_fun(y, x0, P, sz, fun_raw, f0)
%PRECOND_FUN Evaluate a preconditioned objective, gradient, and Hessian.
%
%   [f, g, H] = precond_fun(y, x0, P, sz, fun_raw, f0)
%
%   Inputs
%   ------
%   y : vector
%       Optimization variables in the preconditioned coordinate system.
%
%   x0 : vector
%       Reference point in the original coordinate system.
%
%   P : matrix
%       Preconditioning matrix.
%
%   sz : array
%       Size used to reshape the original optimization variables.
%
%   fun_raw : function handle
%       Function returning the objective, gradient, and Hessian in the
%       original coordinate system.
%
%   f0 : scalar
%       Reference objective value used for scaling.
%
%   Outputs
%   -------
%   f : scalar
%       Scaled objective value.
%
%   g : vector
%       Scaled gradient in the preconditioned coordinates.
%
%   H : matrix
%       Scaled Hessian in the preconditioned coordinates.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Transform optimization variables

    x = x0 + P*y;
    c = reshape(x, sz);


    %% Evaluate objective and derivatives

    [f, gc, Hc] = fun_raw(c);


    %% Scale objective and derivatives

    s = abs(f0);

    f = f / s;
    g = P.' * (gc(:) / s);
    H = P.' * (Hc / s) * P;

    H = 0.5 * (H + H.');

end