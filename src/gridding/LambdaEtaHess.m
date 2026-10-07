function [R, dR, HR] = LambdaEtaHess( ...
    spl, W, D, M, alpha, base, nu_hat, Smax, J, eta)
%LAMBDAETAHESS Compute Lambda, its gradient, and Hessian.
%
%   [R,dR,HR] = LambdaEtaHess( ...
%       spl,W,D,M,alpha,base,nu_hat,Smax,J,eta)
%
%   Inputs
%   ------
%   spl : array
%       Spline coefficients.
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
%   base : char or string
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
%   eta : vector
%       Target error profile.
%
%   Outputs
%   -------
%   R : M x 1 vector
%       Lambda ratio.
%
%   dR : M x L x 2 array
%       Gradient with respect to the real and imaginary coefficients.
%
%   HR : M x 2L x 2L array
%       Hessian with respect to the real and imaginary coefficients.
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
        base = 'spline'
        nu_hat = @(s) double(s==0)
        Smax = 32
        J = []
        eta = []
    end


    %% Input validation

    assert(strcmp(base, 'legendre'), ...
        'Only legendre mode is supported.');

    assert(~isempty(J), ...
        'J must be supplied.');


    %% Dimensions

    [L, ~] = size(spl);
    S = 2*Smax + 1;

    assert( ...
        size(J, 1) == M && size(J, 2) == L && size(J, 3) == S, ...
        'J must be M x L x S.');


    %% Complex spline coefficients

    z = spl(:, 1) + 1i*spl(:, 2);


    %% Fourier coefficients

    Chat = zeros(M, S);

    for s = 1:S
        Chat(:, s) = J(:, :, s) * z;
    end


    %% Frequency range

    s_range = -Smax:Smax;


    %% Weighting vectors and matrices

    nuvec = nu_hat(s_range);
    nuvec = nuvec(:);

    Nnu = nu_hat(s_range - s_range.');
    Nnu = (Nnu + Nnu')/2;

    E = Nnu - nuvec*nuvec';


    %% Outputs

    R  = zeros(M, 1);
    dR = zeros(M, L, 2);
    HR = zeros(M, 2*L, 2*L);


    %% Real coordinate vector

    x = [real(z); imag(z)];


    %% Real block representations

    BN = blkdiag(Nnu, Nnu);
    BE = blkdiag(E, E);


    %% Compute Lambda, gradient, and Hessian

    for m = 1:M

        % -------------------------------------------------------------
        % Construct H so that c = H*z
        % -------------------------------------------------------------

        H = zeros(S, L);

        for s = 1:S
            H(s, :) = J(m, :, s);
        end

        c = Chat(m, :).';


        % -------------------------------------------------------------
        % Realification of complex linear map
        %
        % [Re(c)]       [ Re(H)  -Im(H)] [Re(z)]
        % [Im(c)]   =   [ Im(H)   Re(H)] [Im(z)]
        % -------------------------------------------------------------

        A = real(H);
        B = imag(H);

        T = [A, -B;
             B,  A];


        % -------------------------------------------------------------
        % Denominator
        %
        % v = c' * Nnu * c
        % -------------------------------------------------------------

        v = real(c' * Nnu * c);

        Qv = T.' * BN * T;

        gv = 2*Qv*x;
        Hv = 2*Qv;


        % eta = 0: u = c'*E*c
        % eta = 1: u = -|c'*nu|^2

        a = c' * nuvec;

        r  = T.' * [nuvec; zeros(S, 1)];
        im = T.' * [zeros(S, 1); nuvec];
        K  = r*r.' + im*im.';

        u0  = real(c' * E * c);
        gu0 = 2*(T.'*BE*T)*x;
        Hu0 = 2*(T.'*BE*T);

        u1  = -abs(a)^2;
        gu1 = -2*K*x;
        Hu1 = -2*K;

        u  = (1-eta(m))*u0 + eta(m)*u1;
        gu = (1-eta(m))*gu0 + eta(m)*gu1;
        Hu = (1-eta(m))*Hu0 + eta(m)*Hu1;


        % -------------------------------------------------------------
        % Ratio
        %
        % R = u/v
        % -------------------------------------------------------------

        R(m) = u/v;


        % -------------------------------------------------------------
        % Gradient
        % -------------------------------------------------------------

        d = gu/v - u*gv/v^2;

        dR(m, :, 1) = d(1:L);
        dR(m, :, 2) = d(L+1:end);


        % -------------------------------------------------------------
        % Hessian of u/v
        % -------------------------------------------------------------

        Hm = Hu/v ...
           - (gu*gv' + gv*gu')/v^2 ...
           - u*Hv/v^2 ...
           + 2*u*(gv*gv')/v^3;

        HR(m, :, :) = 0.5*(Hm + Hm.');

    end

end