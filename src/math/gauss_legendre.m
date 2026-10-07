function [nodes, weights] = gauss_legendre(N)
%GAUSS_LEGENDRE Compute N-point Gauss-Legendre quadrature on [-1, 1].
%
%   [nodes, weights] = gauss_legendre(N)
%
%   Inputs
%   ------
%   N : positive integer
%       Number of Gauss-Legendre quadrature points.
%
%   Outputs
%   -------
%   nodes : column vector
%       Gauss-Legendre quadrature nodes on [-1, 1].
%
%   weights : column vector
%       Corresponding Gauss-Legendre quadrature weights.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        N (1,1) double {mustBeInteger, mustBePositive}
    end


    %% Special case: one-point quadrature

    if N == 1

        nodes   = 0;
        weights = 2;

        return

    end


    %% Jacobi matrix

    k = (1:N-1)';

    beta = k ./ sqrt(4*k.^2 - 1);

    % Golub-Welsch: symmetric tridiagonal Jacobi matrix.
    J = diag(beta, 1) + diag(beta, -1);


    %% Eigenvalue decomposition

    [V, D] = eig(J, 'vector');

    [nodes, idx] = sort(D);

    weights = 2 * V(1, idx).^2;


    %% Output orientation

    nodes = nodes(:);
    weights = weights(:);

end