function [f, g, H] = objective_hess(c, W, D, M, alpha, base, ...
    nu_hat, Smax, J, eta, p, rho)
%OBJECTIVE_HESS Compute the objective, gradient, and Hessian.
%
%   [f, g, H] = objective_hess(c, W, D, M, alpha, base, ...
%       nu_hat, Smax, J, eta, p, rho)
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    [R, dR, HR] = LambdaEtaHess( ...
        c, W, D, M, alpha, base, nu_hat, Smax, J, eta);

    JR = [dR(:,:,1), dR(:,:,2)];

    epsR = 0;


    %% Smooth |R|^p

    a = R.^2 + epsR^2;
    phi = a.^(p/2);


    %% Scalarizer

    % -|R|^p + rho*max(R,0)^p
    f = sum(-phi + rho*(max(R,0).^2 + epsR^2).^(p/2));


    %% Derivative

    dfdR = -p*R.*a.^(p/2-1);

    idx = R > 0;
    b = R(idx).^2 + epsR^2;

    dfdR(idx) = dfdR(idx) ...
        + rho*p*R(idx).*b.^(p/2-1);


    %% Second derivative

    d2fdR = -p*a.^(p/2-1) ...
            -p*(p-2)*R.^2.*a.^(p/2-2);

    b = R(idx).^2 + epsR^2;

    d2fdR(idx) = d2fdR(idx) ...
        + rho*p*b.^(p/2-1) ...
        + rho*p*(p-2)*R(idx).^2.*b.^(p/2-2);


    %% Gradient

    g = JR.' * dfdR;
    g = reshape(g, size(c));


    %% Hessian

    H = JR.' * (d2fdR .* JR);

    for m = 1:M
        H = H + dfdR(m) * squeeze(HR(m,:,:));
    end

    H = 0.5 * (H + H.');

end