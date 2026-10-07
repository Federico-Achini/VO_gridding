function c = to_legendre(L, C, N)
%TO_LEGENDRE Compute a Legendre projection of a function on [-1, 1].
%
%   c = to_legendre(L, C, N)
%
%   Inputs
%   ------
%   L : positive integer
%       Number of Legendre coefficients to compute, corresponding to
%       polynomial orders 0 through L-1.
%
%   C : function handle
%       Function C(u) evaluated for u in [-1, 1]. The function may be
%       complex-valued.
%
%   N : positive integer
%       Number of Gauss-Legendre quadrature nodes. Default: 2*L.
%
%   Output
%   ------
%   c : L x 2 array
%       Legendre coefficients. The first column contains the real
%       coefficients and the second column contains the imaginary
%       coefficients.
%
%   Description
%   -----------
%   The coefficients are computed using Gauss-Legendre quadrature
%   projection. This replaces the ill-conditioned equally spaced
%   collocation solve
%
%       M = legendre_poly(L, tau).';
%       c = linsolve(M, fun).
%
%   For a polynomial C of degree <= 2N-1, the quadrature is exact.
%   For a general function, increase N until the coefficients stop
%   changing.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        L (1,1) double {mustBePositive, mustBeInteger}
        C function_handle
        N (1,1) double {mustBePositive, mustBeInteger} = 2*L
    end


    %% Validate quadrature order

    if N < L
        error('legendre_fit:N', 'N must be >= L for a well-posed fit.');
    end


    %% Gauss-Legendre quadrature

    [tau, w] = gauss_legendre(N);      % N x 1 each, well-conditioned nodes
    tau = tau(:).';                    % 1 x N
    w   = w(:).';                      % 1 x N


    %% Legendre basis

    P = legendre_poly(L, tau);         % L x N, P(n+1,:) = P_n(tau)


    %% Function evaluation

    fvals = arrayfun(C, tau);          % 1 x N
    fvals = fvals(:).';


    %% Legendre projection

    n_vec  = (0:L-1).';                % L x 1
    norm_w = (2*n_vec + 1) / 2;        % orthonormality weight per order

    % c_n = (2n+1)/2 * sum_i w_i P_n(tau_i) f(tau_i)
    % (quadrature-weighted projection)
    c_complex = norm_w .* (P * (w .* fvals).');   % L x 1

    c = [real(c_complex), imag(c_complex)];

end