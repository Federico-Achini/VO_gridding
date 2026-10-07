function R = Lambda(spl, W, D, M, alpha, mode, nu_hat, Smax, J)
%LAMBDA Compute the gridding error shape function (Lambda(x)) in 1-D.
%
%   R = Lambda(spl, W, D, M, alpha, mode, nu_hat, Smax, J)
%
%   The function computes the gridding error shape function for an
%   arbitrary time-sample equidistribution measure nu (Hypothesis 1).
%
%   Inputs
%   ------
%   spl : array
%       Kernel representation. Spline structure or sampled kernel
%       depending on the selected mode.
%
%   W : scalar
%       Kernel half-width (support in samples).
%
%   D : scalar
%       Discretization factor (samples per unit).
%
%   M : scalar
%       Number of evaluation points for the output Lambda(x).
%
%   alpha : scalar
%       Oversampling factor gamma. Default: 1.
%
%   mode : char or string
%       Kernel representation mode. Default: 'spline'.
%
%   nu_hat : function handle
%       Fourier coefficients of the limiting time-sample measure nu.
%       The default is @(s) double(s==0), corresponding to uniform
%       time samples.
%
%   Smax : scalar
%       Truncation parameter for the shift sums over s and s'.
%       Default: 32.
%
%   J : array
%       Precomputed Fourier-Legendre basis evaluations.
%
%   Output
%   ------
%   R : vector
%       Values of the relative gridding error.
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
    end


    %% Time-domain sample grid and sampled kernel

    if strcmp(mode, 'spline')
        dd = 100;
        dt = 1/(D*dd);
        N = 2*W*D;

        tgrid = discretize(-W, W, dd*N, 'both');
        ker = eval_spline_fast(spl, tgrid);

        tgrid = orientation(tgrid, 'h');
        ker   = orientation(ker,   'h');

    elseif strcmp(mode, 'vector')
        dt = 1/D;
        N = 2*W*D;

        tgrid = discretize(-W, W, N, 'both');
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


    %% Gridding-error component

    E = Nnu - nuvec * nuvec';


    %% Error evaluation

    num = sum((Chat * E) .* conj(Chat), 2);
    den = sum((Chat * Nnu) .* conj(Chat), 2);

    R = real(num ./ den);


    %% Remove insignificant negative values caused by roundoff

    tol = 100 * eps(max(1, max(abs(R))));
    R(R < 0 & R > -tol) = 0;

end