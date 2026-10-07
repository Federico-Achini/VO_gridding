function R = LambdaEta(spl, W, D, M, alpha, mode, nu_hat, Smax, J, eta)
%LAMBDAETA Compute the gridding error shape function Lambda(x).
%
%   R = LambdaEta(spl, W, D, M, alpha, mode, nu_hat, Smax, J, eta)
%
%   Computes the gridding error shape function Lambda(x) for a given
%   kernel and time-sample equidistribution measure nu.
%
%   Inputs
%   ------
%   spl : array
%       Kernel representation.
%
%   W : scalar
%       Kernel half-width.
%
%   D : scalar
%       Discretization factor.
%
%   M : scalar
%       Number of evaluation points for the output Lambda(x).
%
%   alpha : scalar
%       Oversampling factor. Default: 1.
%
%   mode : char or string
%       Kernel representation mode. Default: 'spline'.
%
%   nu_hat : function handle
%       Fourier coefficients of the limiting time-sample measure nu.
%       Default: @(s) double(s==0).
%
%   Smax : scalar
%       Maximum shift used in the truncated Fourier sums.
%       Default: 32.
%
%   J : array
%       Precomputed Fourier-Legendre basis evaluations.
%
%   eta : vector
%       Target error profile.
%
%   Output
%   ------
%   R : vector
%       Values of the relative gridding error Lambda(x).
%
%   ----------------------------------------------------------------------
%
%   The implemented expression is
%
%       ell*(w; C) = 1 - |sum_s nu_hat(s) *
%                         conj(Chat(w/gamma - s))|^2
%                    / sum_{s,s'} nu_hat(s'-s) *
%                         Chat(w/gamma-s) *
%                         conj(Chat(w/gamma-s')).
%
%   The (1+o(1)) factor from the asymptotic bound is omitted.
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
        Smax = 32
        J = []
        eta = []
    end


    %% Time-domain sample grid and sampled kernel

    if strcmp(mode, 'spline')
        dd = 51;
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
        %

    else
        error('Insert mode is not allowed.');
    end


    %% x-domain for evaluation

    % M points in [-0.5, 0.5].
    x = discretize(-0.5, 0.5, M);
    x = orientation(x, 'v');


    %% Fourier coefficients of the kernel

    % Chat(:,k) = Chat(x/alpha - s_range(k)) for every x.
    s_range = -Smax:Smax;

    if strcmp(mode, 'legendre')

        if isempty(J)
            S = 2*Smax + 1;
            Chat = zeros(M, S);

            for s = 1:S
                J = 2*W*j_basis(L, W*(x/alpha - s_range(s)));
                Chat(:, s) = 2*W*J.' * ...
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


    %% Time-sample measure

    nuvec = nu_hat(s_range);
    nuvec = nuvec(:);

    Nnu = nu_hat(s_range - s_range.');

    % Symmetrize to suppress tiny numerical non-Hermiticity.
    Nnu = (Nnu + Nnu')/2;


    %% Error evaluation

    % Evaluate the error and total quadratic forms directly.
    den = sum((Chat * Nnu) .* conj(Chat), 2);

    idx0 = (eta == 0);
    idx1 = (eta == 1);

    R = zeros(M, 1);


    %% eta = 0

    E = Nnu - nuvec * nuvec';

    num0 = sum( ...
        (Chat(idx0, :) * E) .* conj(Chat(idx0, :)), ...
        2);

    R(idx0) = real(num0 ./ den(idx0));


    %% eta = 1

    a1 = Chat(idx1, :) * nuvec;

    R(idx1) = -abs(a1).^2 ./ den(idx1);

end