%% 
% ---------------
% |  EXAMPLE.m  |
% ---------------
%
% This script optimizes and evaluates a gridding interpolation kernel
% for different frequency-distribution widths.
%
% The workflow includes:
%   - MATLAB environment and path initialization.
%   - Definition of kernel and optimization parameters.
%   - Construction of a Kaiser-Bessel initial kernel.
%   - Computation of the initial gridding error functional.
%   - Kernel optimization for multiple frequency-distribution widths.
%   - Evaluation of the optimized kernel and its error functional.
%   - Comparison with the Kaiser-Bessel reference kernel.
%   - Reconstruction tests using synthetic signals.
%   - Measurement of gridding computation times.
%   - Analysis and visualization of reconstruction errors.
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

p = 0.5;

random_signals = true;
use_real_kernel = true;


%% Target error profile

x = discretize(-0.5,0.5,M);
eta = zeros(size(x));

spl0 = spl_kaiser(W,D);

shift = @(p,x) exp(2i*pi*p*x);
guess = @(x) eval_spline_fast(spl0,W*x).*shift(0.00*W,x);
c00 = to_legendre(L,guess,2*W*D+1);

nu_hat_0 = @(s) double (s == 0);
Smax = 64;

J = create_J_matrix(L,W,x,alpha,-Smax:Smax);

h0 = correction(c00,W,D,M,alpha,base,nu_hat_0,Smax,J);
sp0 = prepare_spline_fast(spl0, W, D);
Lambda_initial = Lambda(c00, W, D, M, alpha, base, nu_hat_0, Smax, J);


%% Frequency-distribution parameter sweep

s_nus = [1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 5e-1, 1, 5, 10];
s_nus = s_nus(1:end);
s_l = length(s_nus);

ctime_vo = zeros(size(s_nus));
ctime_kb = zeros(size(s_nus));

err_vo = zeros(s_l, M);
err_kb = zeros(s_l, M);


for ss = 1:s_l

    s_nu = s_nus(ss);
    nu_hat = @(s) exp(-2*pi^2*s_nu^2*alpha^2*s.^2);
    Smax = max(ceil(3 / (2 * alpha * pi * s_nu)), 64);

    fprintf("Test number %d\n - sigma = %.3f\n - Smax = %d\n", ...
        ss, s_nu, Smax);

    J = create_J_matrix(L,W,x,alpha,-Smax:Smax);

    [R_initial, ~, ~] = LambdaEtaHess( ...
        c00,W,D,M,alpha,base,nu_hat,Smax,J,eta);

    A = norm(R_initial,p)^p;
    B = norm(max(R_initial,0),p)^p;

    B_floor = max(1e-12 * A, B);
    rho = 1.5*A/B_floor;

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


    %% Plot optimized and reference kernels

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


    %% Define testing parameters

    Ns = [1e3, 1e7];
    Ns_size = length(Ns);
    trials = 30;


    %% Test the optimized kernel

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


        %% Generate test signals

        if random_signals

            signals = randn(N, trials); % + 1i*randn(N, trials);
            signals = signals ./ vecnorm(signals, 2, 1);
            phases = mod((0:N-1).' + randn(N, trials)*s_nu, N);

        end


        %% Run reconstruction trials

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


        %% Store results

        err_vo(ss, :) = mean(tmp_err_vo, 1);
        err_kb(ss, :) = mean(tmp_err_kb, 1);

        ctime_vo(ss) = mean(tmp_time_vo);
        ctime_kb(ss) = mean(tmp_time_kb);

    end

end


%% Plot reconstruction errors

f_nus = figure('Name','Nus');

semilogy(x, err_kb(1, :), ...
    'Color',blue(1), ...
    'DisplayName','KB', ...
    'LineWidth',2.5, ...
    'LineStyle',':');

hold on

for ss = s_l-3:-1:1

    semilogy(x, err_vo(ss, :), ...
        'DisplayName',strcat( ...
            'VO ($\sigma = ', string(s_nus(ss)), '$)'), ...
        'LineWidth',2.5);

end

legend('Location','best')


%% Compute mean reconstruction errors

mean_errs_vo = zeros(s_l - 3, 1);

for ss = 1:s_l - 3
    mean_errs_vo(ss) = mean(err_vo(ss, :));
end

mean_errs_kb = mean(err_kb(ss, :));

