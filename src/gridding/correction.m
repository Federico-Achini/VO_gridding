function y = correction(spl, W, D, M, alpha, mode, nu_hat, Smax, J)
%CORRECTION Compute the correction function for the gridding kernel.
%
%   y = correction(spl, W, D, M, alpha, mode, nu_hat, Smax, J)
%
%   Inputs
%   ------
%   spl : array or spline representation
%       Kernel representation.
%
%   W : scalar
%       Kernel width.
%
%   D : scalar
%       Kernel discretization parameter.
%
%   M : scalar
%       Number of evaluation points.
%
%   alpha : scalar
%       Oversampling factor.
%
%   mode : char or string
%       Kernel representation mode.
%
%   nu_hat : function handle
%       Fourier-domain weighting function.
%
%   Smax : scalar
%       Maximum Fourier shift.
%
%   J : array
%       Precomputed Fourier-Legendre basis evaluations.
%
%   Output
%   ------
%   y : vector
%       Correction function evaluated at the M points.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        spl
        W
        D
        M
        alpha = 1
        mode = 'spline'
        nu_hat = @(s) double(s==0)
        Smax = 1
        J = []
    end


    %% Time-domain sample grid and sampled kernel

    if strcmp(mode, 'spline')
        dd = 115;
        dt = 1/(D*dd);
        N = 2*W*D;

        tgrid = discretize(-W, W, dd*N + 1, 'both');
        ker = eval_spline_fast(spl, tgrid);

        tgrid = orientation(tgrid, 'h');
        ker   = orientation(ker,   'h');

    elseif strcmp(mode, 'vector')
        N = 2*W*D;

        tgrid = discretize(-W, W, N + 1, 'both');
        ker = spl;

        tgrid = orientation(tgrid, 'h');
        ker   = orientation(ker,   'h');

    elseif strcmp(mode, 'legendre')
        [L, ~] = size(spl);

    else
        error('Insert mode is not allowed.');
    end


    %% x-domain for evaluation

    % M points in [-0.5, 0.5].
    x = discretize(-0.5, 0.5, M);
    x = orientation(x, 'v');


    %% Fourier-Legendre basis evaluation

    % Chat(x/alpha - s) for every evaluation point x and every shift s.
    s_range = -Smax:Smax;

    if strcmp(mode, 'legendre')

        if isempty(J)
            S = 2*Smax + 1;
            Chat = zeros(M, S);

            for s = 1:S
                Jtmp = j_basis(L, W * (x/alpha - s_range(s)));
                Chat(:, s) = 2*W * Jtmp.' * ...
                    (spl(:, 1) + 1i*spl(:, 2));
            end

        else
            [~, ~, S] = size(J);

            assert( ...
                S == 2*Smax + 1, ...
                'The values of the Fourier truncation must be consistent.');

            Chat = zeros(M, S);

            for s = 1:S
                Chat(:, s) = J(:, :, s) * ...
                    (spl(:, 1) + 1i*spl(:, 2));
            end
        end

    else
        Fx = exp(-2i*pi * (x/alpha) * tgrid);
        Gs = exp(2i*pi * tgrid.' * s_range);

        Chat = dt * (Fx .* ker) * Gs;
    end


    %% Weighting vectors and matrix

    nuvec = nu_hat(s_range);
    nuvec = nuvec(:);

    Nnu = nu_hat(s_range - s_range.');

    % Enforce Hermitian structure numerically.
    Nnu = (Nnu + Nnu')/2;


    %% Correction function

    num = conj(Chat) * nuvec;

    den = sum((Chat * Nnu) .* conj(Chat), 2);
    den = real(den);

    y = num ./ den;

end