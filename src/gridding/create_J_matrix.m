function J = create_J_matrix(L, W, x, alpha, s_range)
%CREATE_J_MATRIX Precompute shifted Fourier-Legendre basis evaluations.
%
%   J = create_J_matrix(L, W, x, alpha, s_range)
%
%   Inputs
%   ------
%   L : scalar
%       Number of Fourier-Legendre basis functions.
%
%   W : scalar
%       Kernel width.
%
%   x : vector
%       Evaluation points.
%
%   alpha : scalar
%       Scaling factor.
%
%   s_range : vector
%       Range of Fourier shifts.
%
%   Output
%   ------
%   J : array
%       Precomputed shifted Fourier-Legendre basis evaluations.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        L
        W
        x
        alpha
        s_range
    end


    %% Dimensions

    N = length(x);
    S = length(s_range);

    J = zeros(L, N, S);


    %% Fourier-Legendre basis evaluations

    for i = 1:S
        J(:, :, i) = ...
            2*W*j_basis(L, W*(x/alpha - s_range(i)));
    end


    %% Rearrange dimensions

    J = permute(J, [2 1 3]);

end