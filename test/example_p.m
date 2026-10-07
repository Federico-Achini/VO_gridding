% --------------------
% |  EXAMPLE.m       |
% --------------------
%
% This script investigates the effect of the scalarization exponent p
% on kernel optimization for gridding error control.
%
% The workflow defines a spatially localized target error profile,
% optimizes the gridding kernel for several values of p, and evaluates
% the resulting kernels using reconstruction errors from random signals.
% The optimized kernels are compared against a Kaiser-Bessel reference.
%
% Author: Federico Achini
% Date:   August 2025
%
% ------------------------------------------------------------------------

%% Add necessary paths and setup

clc;
startup
addpath('C:\Users\FedericoAchini\OneDrive - Revenio Group Oyj\Matlab\finufft\finufft\build\matlab')
%mip load finufft


%% Parameters

W = 3;
D = 300;
M = 100;
L = 15*W;
alpha = 2;

compute_kernel = true;
base = 'legendre';
compute_times = true;

random_signals = true;
use_real_kernel = false;


%% Target error profile

x = discretize(-0.5,0.5,M);
eta = ones(size(x));
eta(abs(x - 0.25) < 0.05) = 0;

spl0 = spl_kaiser(W,D);

shift = @(p,x) exp(2i*pi*p*x);
guess = @(x) eval_spline_fast(spl0,W*x).*shift(0.00*W,x);
c00 = to_legendre(L,guess,2*W*D+1);

nu_hat_0 = @(s) double (s == 0);
nu_hat = nu_hat_0;
Smax = 64;

J = create_J_matrix(L,W,x,alpha,-Smax:Smax);

h0 = correction(c00,W,D,M,alpha,base,nu_hat_0,Smax,J);
sp0 = prepare_spline_fast(spl0, W, D);

Lambda_initial = Lambda( ...
    c00, W, D, M, alpha, base, nu_hat_0, Smax, J);

[R_initial, ~, ~] = LambdaEtaHess( ...
    c00,W,D,M,alpha,base,nu_hat,Smax,J,eta);


%% Scalarization parameter sweep

ps = [0.1, 0.25, 0.5, 0.75, 1, 1.5, 2];
pl = length(ps);

ctime_vo = zeros(size(ps));
ctime_kb = zeros(size(ps));

err_vo = zeros(pl, M);
err_kb = zeros(pl, M);


for pp = 1:pl

    p = ps(pp);

    fprintf("Test number %d\n - p = %.3f\n", pp, p);

    rho = 1e1;

    s0 = abs(scalarizer(R_initial,p,rho));

    c_current = c00;
    nIter = 3;


    %% Kernel optimization

    for i = 1:nIter

        c0 = c_current;
        n  = numel(c0);

        fun_raw = @(c) objective_hess( ...
            c,W,D,M,alpha,base,nu_hat,Smax,J,eta,p,rho);

        [f0,~,H0] = fun_raw(c0);

        x0 = c0(:);
        e  = x0/norm(x0);
        Z = null(e.');

        H0r = Z.'*H0*Z;
        H0r = (H0r + H0r.')/2;

        [V,Dd] = eig(H0r);
        lam = diag(Dd);
        lamMax = max(abs(lam));

        P = Z*V*diag(sqrt(abs(f0))./sqrt(abs(lam)));

        fun = @(y) precond_fun( ...
            y,x0,P,size(c0),fun_raw,f0);


        opts = optimoptions('fminunc', ...
            'Algorithm','trust-region', ...
            'SpecifyObjectiveGradient',true, ...
            'HessianFcn','objective', ...
            'Display','iter', ...
            'MaxIterations',7e1, ...
            'MaxFunctionEvaluations',1e6, ...
            'OptimalityTolerance',0, ...
            'StepTolerance',1e-6, ...
            'FunctionTolerance',1e-16);

        [y,fval,exitflag,output] = fminunc( ...
            fun,zeros(n-1, 1),opts);


        % Recover c
        x_min = x0 + P*y;
        c_min = reshape(x_min,size(c0));
        c_current = c_min;

    end


    %% Normalize optimized kernel

    c_min = c_min ./ norm(c_min);

    if use_real_kernel
        c_min(:, 2) = 0;
    end

    spl = create_legendre_spline(c_min, W, D);


    %% Evaluate error functional

    Lambda_final = Lambda( ...
        c_min, W, D, M, alpha, base, nu_hat, Smax, J);

    f_lambda = figure('Name', 'Lambda');

    semilogy(x, Lambda_final, ':x', ...
        'DisplayName', '$\Lambda_{fin}$', ...
        'Color', red(0))

    hold on

    semilogy(x, Lambda_initial, ':o', ...
        'DisplayName', '$\Lambda_{in}$', ...
        'Color', blue(0))

    title("Error functional $\Lambda$")
    legend('Location','best');
    xlabel("$\omega$")
    ylabel("$\Lambda(\omega)$")


    %% Plot kernel

    tau = discretize(-W, W, 2*W*D*10 + 1, 'both');

    ker = eval_spline_fast(spl, tau);
    ker0 = eval_spline_fast(spl0, tau);

    f_kerr = figure('Name', 'Kernel real');

    plot(tau, real(ker), 'Color', blue(0))
    hold on
    plot(tau, real(ker0), 'Color', blue(1))

    f_keri = figure('Name', 'Kernel real');

    plot(tau, imag(ker), 'Color', red(0))
    hold on
    plot(tau, imag(ker0), 'Color',red(1))


    %% Compute gridding parameters

    h = correction(c_min,W,D,M,alpha,base,nu_hat,Smax,J);
    sp = prepare_spline_fast(spl, W, D);


    %% Testing parameters

    Ns = [1e3]; %, 1e7];
    Ns_size = length(Ns);
    trials = 30;


    %% Testing the optimized kernel

    for d = 1:Ns_size

        N = Ns(d);

        if N > 1e5

            compute_true = false;
            tmp_time_vo = zeros(trials, 1);
            tmp_time_kb = zeros(trials, 1);

        else

            compute_true = true;
            tmp_err_vo = zeros(trials, M);
            tmp_err_kb = zeros(trials, M);

        end


        %% Generate random signals

        if random_signals

            signals = randn(N, trials); % + 1i*randn(N, trials);
            signals = signals ./ vecnorm(signals, 2, 1);
            phases = mod((0:N-1).' + randn(N, trials)*N, N);

        end


        %% Run multiple trials

        for i = 1:trials

            fprintf('Trial number %d\n', i);

            % Generate signal based on test mode
            if random_signals

                signal = signals(:, i);
                phase = phases(:, i);

            else

                phase = mod((0:N-1).' + randn(N, 1)*s_nu, N);
                R = randi(1) + 1;
                refl = rand(R, 1) - 0.5;
                ampl = abs(randn(R, 1));
                signal = gen_signal_from_refl(phase, refl, ampl, true);

            end


            if compute_true

                % True reconstruction (ground truth)
                y_true = reconstruction(phase, signal, x);

            end


            tic

            y_kb = gridding_al(signal, phase, M, W, sp0, h0, alpha);

            tmp_time_kb(i) = toc;


            tic

            y_vo = gridding_al(signal, phase, M, W, sp, h, alpha);

            tmp_time_vo(i) = toc;


            if compute_true

                tmp_err_vo(i, :) = abs(y_true - y_vo);
                tmp_err_kb(i, :) = abs(y_true - y_kb);

            end

        end


        err_vo(pp, :) = mean(tmp_err_vo, 1);
        err_kb(pp, :) = mean(tmp_err_kb, 1);

        ctime_vo(pp) = mean(tmp_time_vo);
        ctime_kb(pp) = mean(tmp_time_kb);

    end

end


%% Plot reconstruction errors

f_nus = figure('Name','Nus');

semilogy(x, err_kb(1, :), ...
    'Color',blue(1), 'DisplayName','KB', 'LineWidth',2.5, 'LineStyle',':');

hold on

for pp = 1:pl

    semilogy(x, err_vo(pp, :), ...
        'DisplayName',strcat('VO ($p = ', string(ps(pp)), '$)'), ...
        'LineWidth',2.5);

end

legend('Location','best')


%% Compute mean reconstruction errors

mean_errs_vo = zeros(pl, 1);

for pp = 1:pl
    mean_errs_vo(pp) = mean(err_vo(pp, :));
end

mean_errs_kb = mean(err_kb(pp, :));

