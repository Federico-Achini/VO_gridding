%% ------------------------------------------------------------------------
%  EXAMPLE.m
%  ------------------------------------------------------------------------
%
%  This script demonstrates kernel optimization for gridding error control.
%
%  Author: Federico Achini
%  Date:   August 2025
%
%  ------------------------------------------------------------------------

%% Add necessary paths and setup

clc;
startup
mip load finufft


%% ---- Parameters --------------------------------------------------------

W = 3;
D = 300;
M = 100;
L = 15 * W;
alpha = 2;

compute_kernel   = true;
base             = 'legendre';
test_name        = 3;
compute_analysis = false;
compute_times    = true;

p = 0.5;

random_signals = true;
s_nu = 10;

nu_hat = @(s) exp(-2 * pi^2 * s_nu^2 * alpha^2 * s.^2);
Smax = 64;

use_real_kernel = false;

x = discretize(-0.5, 0.5, M);

J = create_J_matrix( ...
    L, W, x, alpha, -Smax:Smax);


%% ---- Target error profile ----------------------------------------------

eta = ones(size(x));

switch test_name

    case 0

        eta = zeros(size(eta));

    case 1

        eta(x > 0) = 0;

    case 2

        eta(abs(x - 0.25) < 0.05) = 0;

    case 3

        eta(abs(x + 0.3) < 0.05) = 1e-25;
        eta(abs(x - 0.1) < 0.1)  = 1e-27;

end


f_eta = figure('Name', 'eta');

plot(x, eta);

title("");
xlabel("$\omega$");
ylabel("$\eta(\omega)$");
ylim([-1e-1, 1.1]);


%% ---- Initial guess -----------------------------------------------------

spl0 = spl_kaiser(W, D);

shift = @(p, x) exp(2i * pi * p * x);
guess = @(x) eval_spline_fast(spl0, W * x) .* shift(0.00 * W, x);

c00 = to_legendre(L, guess, 2 * W * D + 1);

spl0_comp = create_legendre_spline(c00, W, D);


%% ---- Initial Lambda ----------------------------------------------------

Lmb_opt_coef = Lambda( ...
    c00, W, D, M, alpha, base, nu_hat, Smax, J);

Lmb_opt_spl_comp = Lambda( ...
    spl0_comp, W, D, M, alpha, 'spline', nu_hat, Smax, J);

Lmb_opt_spl_orig = Lambda( ...
    spl0, W, D, M, alpha, 'spline', nu_hat, Smax, J);


figure;

semilogy( ...
    x, Lmb_opt_coef, ':xr', ...
    'DisplayName', 'leg');

hold on;

semilogy( ...
    x, Lmb_opt_spl_comp, ':*g', ...
    'DisplayName', 'spl comp');

semilogy( ...
    x, Lmb_opt_spl_orig, '-ob', ...
    'DisplayName', 'spl orig');

legend;
grid on;

xlabel('x');
ylabel('$\Lambda(x)$');


%% ---- Initial residual --------------------------------------------------

[R_initial, ~, ~] = LambdaEtaHess( ...
    c00, W, D, M, alpha, base, nu_hat, Smax, J, eta);


%% ---- Set rho -----------------------------------------------------------

A = norm(R_initial, p)^p;
B = norm(max(R_initial, 0), p)^p;

B_floor = max(1e-12 * A, B);
rho = 1e8 * A / B_floor;

fprintf( ...
    'A = %.6e, B = %.6e, rho*B = %.6e\n', ...
    A, B, rho * B);

s0 = abs(scalarizer(R_initial, p, rho));

if ~isfinite(s0) || s0 == 0
    error('Invalid initial scalarization.');
end

fprintf('rho = %.6e\n', rho);


%% ---- fminunc + Hessian eigenvalue preconditioning ---------------------

[R, dR, HR] = LambdaEtaHess( ...
    c00, W, D, M, alpha, base, nu_hat, Smax, J, eta);

JR = [dR(:, :, 1), dR(:, :, 2)];

g_A = JR.' * (-2 * R);

idx = R > 0;

dB_dR = zeros(size(R));
dB_dR(idx) = 2 * R(idx);

g_B = JR.' * dB_dR;

fprintf('||g_A|| = %.16e\n', norm(g_A));
fprintf('||g_B|| = %.16e\n', norm(g_B));
fprintf('||g_A+rho*g_B|| = %.16e\n', ...
    norm(g_A + rho * g_B));


c_current = c00;
nIter = 6;

for i = 1:nIter

    c0 = c_current;
    n = numel(c0);

    fun_raw = @(c) objective_hess( ...
        c, W, D, M, alpha, base, nu_hat, Smax, J, eta, p, rho);

    [f0, ~, H0] = fun_raw(c0);


    % ---------------------------------------------------------------------
    % Remove homogeneous direction
    % ---------------------------------------------------------------------

    x0 = c0(:);
    e = x0 / norm(x0);

    Z = null(e.');


    % ---------------------------------------------------------------------
    % Hessian in scale-free coordinates
    % ---------------------------------------------------------------------

    H0r = Z.' * H0 * Z;
    H0r = (H0r + H0r.') / 2;


    % ---------------------------------------------------------------------
    % Hessian eigenvalue preconditioning
    % ---------------------------------------------------------------------

    [V, Dd] = eig(H0r);
    lam = diag(Dd);

    fprintf('min |lambda| = %.6e\n', min(abs(lam)));
    fprintf('max |lambda| = %.6e\n', max(abs(lam)));
    fprintf( ...
        'cond         = %.6e\n', ...
        max(abs(lam)) / min(abs(lam)));

    lamMax = max(abs(lam));

    P = Z * V * ...
        diag(sqrt(abs(f0)) ./ sqrt(abs(lam)));

    [~, g0, H0check] = fun_raw(c0);

    H0y = P.' * (H0check / abs(f0)) * P;
    H0y = (H0y + H0y.') / 2;

    eig_y = eig(H0y);

    fun = @(y) precond_fun( ...
        y, x0, P, size(c0), fun_raw, f0);


    opts = optimoptions( ...
        'fminunc', ...
        'Algorithm', 'trust-region', ...
        'SpecifyObjectiveGradient', true, ...
        'HessianFcn', 'objective', ...
        'Display', 'iter-detailed', ...
        'MaxIterations', 2e2, ...
        'MaxFunctionEvaluations', 1e6, ...
        'OptimalityTolerance', 0, ...
        'StepTolerance', 1e-8, ...
        'FunctionTolerance', 1e-16);


    [y, fval, exitflag, output] = fminunc( ...
        fun, zeros(n - 1, 1), opts);


    % Recover c

    x_min = x0 + P * y;
    c_min = reshape(x_min, size(c0));

    c_current = c_min;

end


% Evaluate the physical solution

[f_final, g_final, H_final] = fun_raw(c_min);

fprintf('f = %.16e\n', f_final);
fprintf('iterations = %d\n', output.iterations);
fprintf('exitflag = %d\n', exitflag);


%% ---- Final -------------------------------------------------------------

fprintf('Final objective = %.16g\n', f_final);
fprintf( ...
    'Final gradient norm = %.6e\n', ...
    norm(g_final(:)));


c_min = c_min ./ norm(c_min);

if use_real_kernel
    c_min(:, 2) = 0;
end

[R_final, ~, ~] = LambdaEtaHess( ...
    c_min, W, D, M, alpha, base, nu_hat, Smax, J, eta);


%% ---- Optimal kernel evaluation -----------------------------------------

L_test = L;
W_test = W + 1;
D_test = D;
Smaxm = 4 * Smax;

Jm = create_J_matrix( ...
    L_test, W_test, x, alpha, -Smaxm:Smaxm);

spl_test = spl_kaiser(W_test, D_test);

c_test = spline_to_legendre( ...
    spl_test, L_test, W_test);

spl_test_comp = create_legendre_spline( ...
    c_test, W_test, D_test);


%% ---- Lambdas -----------------------------------------------------------

Lambda_initial = Lambda( ...
    c00, W, D, M, alpha, base, nu_hat, Smax, J);

Lambda_test = Lambda( ...
    c_test, W_test, D_test, M, alpha, base, nu_hat, Smaxm, Jm);

Lambda_final = Lambda( ...
    c_min, W, D, M, alpha, base, nu_hat, Smax, J);


% Compute spline Lambdas

spl = create_legendre_spline(c_min, W, D);

Lambda_initial_spline = Lambda( ...
    spl0_comp, W, D, M, alpha, 'spline', nu_hat, Smax);

Lambda_test_spline = Lambda( ...
    spl_test, W_test, D_test, M, alpha, 'spline', nu_hat, Smaxm);

Lambda_final_spline = Lambda( ...
    spl, W, D, M, alpha, 'spline', nu_hat, Smax);


%% ---- Plot Lambdas ------------------------------------------------------

f_lambda = figure('Name', 'Lambda');

semilogy( ...
    x, Lambda_final, ':x', ...
    'DisplayName', '$\Lambda_{fin}$', ...
    'Color', red(0));

hold on;

semilogy( ...
    x, Lambda_initial, ':o', ...
    'DisplayName', '$\Lambda_{in}$', ...
    'Color', blue(0));

semilogy( ...
    x, Lambda_test, ':o', ...
    'DisplayName', '$\Lambda_{opt}$', ...
    'Color', green(0));

semilogy( ...
    x, Lambda_final_spline, '-', ...
    'DisplayName', '$\Lambda_{fin}$', ...
    'Color', red(1));

semilogy( ...
    x, Lambda_initial_spline, '-', ...
    'DisplayName', '$\Lambda_{in}$', ...
    'Color', blue(1));

semilogy( ...
    x, Lambda_test_spline, '-', ...
    'DisplayName', '$\Lambda_{opt}$', ...
    'Color', green(1));

semilogy( ...
    x, eta, ':', ...
    'DisplayName', '$\eta$', ...
    'Color', [.7 .7 .7]);

title("Error functional $\Lambda$");

legend('Location', 'best');

xlabel("$\omega$");
ylabel("$\Lambda(\omega)$");


%% ---- Analysis ----------------------------------------------------------

if compute_analysis

    % Basic measures

    metrics = struct();

    metrics.initial_max_abs = max(abs(R_initial));
    metrics.final_max_abs = max(abs(R_final));

    metrics.initial_L1 = mean(abs(R_initial));
    metrics.final_L1 = mean(abs(R_final));

    metrics.initial_L2 = norm(R_initial, 2) / sqrt(M);
    metrics.final_L2 = norm(R_final, 2) / sqrt(M);

    metrics.initial_Linf = norm(R_initial, inf);
    metrics.final_Linf = norm(R_final, inf);

    metrics.initial_max_pos = max(R_initial);
    metrics.final_max_pos = max(R_final);

    metrics.initial_max_neg = min(R_initial);
    metrics.final_max_neg = min(R_final);

    metrics.initial_pos_fraction = mean(R_initial > 0);
    metrics.final_pos_fraction = mean(R_final > 0);

    metrics.initial_integral_abs = ...
        trapz(x, abs(R_initial));

    metrics.final_integral_abs = ...
        trapz(x, abs(R_final));


    fprintf('\n----- Error metrics -----\n');
    fprintf('                 Initial          Final\n');

    fprintf( ...
        'L1 mean       %12.4e   %12.4e\n', ...
        metrics.initial_L1, ...
        metrics.final_L1);

    fprintf( ...
        'L2 RMS        %12.4e   %12.4e\n', ...
        metrics.initial_L2, ...
        metrics.final_L2);

    fprintf( ...
        'Linf          %12.4e   %12.4e\n', ...
        metrics.initial_Linf, ...
        metrics.final_Linf);

    fprintf( ...
        'Max positive  %12.4e   %12.4e\n', ...
        metrics.initial_max_pos, ...
        metrics.final_max_pos);

    fprintf( ...
        'Min negative  %12.4e   %12.4e\n', ...
        metrics.initial_max_neg, ...
        metrics.final_max_neg);

    fprintf( ...
        'Positive frac %12.4f   %12.4f\n', ...
        metrics.initial_pos_fraction, ...
        metrics.final_pos_fraction);

    fprintf( ...
        'Integral |R|  %12.4e   %12.4e\n', ...
        metrics.initial_integral_abs, ...
        metrics.final_integral_abs);


    %% ---- Error improvement ---------------------------------------------

    improvement = ...
        abs(R_initial) ./ max(abs(R_final), eps);

    figure;

    semilogy( ...
        x, improvement, '-', ...
        'DisplayName', 'Error reduction factor');

    grid on;
    legend('Location', 'best');

    xlabel('$x$', 'Interpreter', 'latex');

    ylabel( ...
        '$|R_{\rm initial}|/|R_{\rm final}|$', ...
        'Interpreter', 'latex');

    title('Pointwise error reduction');


    %% ---- Positive error ------------------------------------------------

    figure;

    semilogy( ...
        x, max(R_initial, 0), '--', ...
        'DisplayName', 'Initial');

    hold on;

    semilogy( ...
        x, max(R_final, 0), '-', ...
        'DisplayName', 'Optimized');

    grid on;
    legend('Location', 'best');

    xlabel('$x$', 'Interpreter', 'latex');

    ylabel( ...
        '$\max(\Lambda-\eta,0)$', ...
        'Interpreter', 'latex');

    title('Positive error component');


    %% ---- Negative error ------------------------------------------------

    figure;

    semilogy( ...
        x, max(-R_initial, 0), '--', ...
        'DisplayName', 'Initial');

    hold on;

    semilogy( ...
        x, max(-R_final, 0), '-', ...
        'DisplayName', 'Optimized');

    grid on;
    legend('Location', 'best');

    xlabel('$x$', 'Interpreter', 'latex');

    ylabel( ...
        '$\max(\eta-\Lambda,0)$', ...
        'Interpreter', 'latex');

    title('Negative error component');


    %% ---- Coefficient change --------------------------------------------

    dc = c_min - c00;

    figure;

    semilogy( ...
        1:L, abs(dc(:, 1)), '-', ...
        'DisplayName', 'Real coefficients');

    hold on;

    semilogy( ...
        1:L, abs(dc(:, 2)), '--', ...
        'DisplayName', 'Imaginary coefficients');

    grid on;
    legend('Location', 'best');

    xlabel('Basis index');

    ylabel( ...
        '$|\Delta c|$', ...
        'Interpreter', 'latex');

    title('Change in optimization coefficients');

end


%% ---- Plot kernel -------------------------------------------------------

tau = discretize( ...
    -W, W, 2 * W * D * 10 + 1, 'both');

ker = eval_spline_fast(spl, tau);
ker0 = eval_spline_fast(spl0, tau);


f_kerr = figure('Name', 'Kernel real');

plot( ...
    tau, real(ker), ...
    'Color', blue(0));

hold on;

plot( ...
    tau, real(ker0), ...
    'Color', blue(1));

hold off;


f_keri = figure('Name', 'Kernel real');

plot( ...
    tau, imag(ker), ...
    'Color', red(0));

hold on;

plot( ...
    tau, imag(ker0), ...
    'Color', red(1));

hold off;


%% ---- Compute correcting functions -------------------------------------

h = correction( ...
    c_min, W, D, M, alpha, base, nu_hat, Smax, J);

h0 = correction( ...
    c00, W, D, M, alpha, base, nu_hat, Smax, J);

h_test = correction( ...
    c_test, W_test, D_test, M, alpha, base, nu_hat, Smaxm, Jm);


%% ---- Spline preparation for algorithm ---------------------------------

sp = prepare_spline_fast( ...
    spl, W, D);

sp0 = prepare_spline_fast( ...
    spl0, W, D);

sp_test = prepare_spline_fast( ...
    spl_test, W_test, D_test);


%% ------------------------------------------------------------------------

Ns = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8];
Ns_size = length(Ns);

ctime_vo = zeros(Ns_size, 1);
ctime_kb = zeros(Ns_size, 1);
ctime_fi = zeros(Ns_size, 1);


%% ---- Testing the optimized kernel -------------------------------------

for d = 1:Ns_size

    N = Ns(d)

    if N > 1e5
        compute_true = false;
    else
        compute_true = true;
    end

    trials = 10;


    % Initialize error storage

    err = zeros(trials, M);
    err_test = zeros(trials, M);
    err_fi = zeros(trials, M);
    err_0 = zeros(trials, M);

    time = zeros(trials, 1);
    time_test = zeros(trials, 1);
    time_fi = zeros(trials, 1);
    time_0 = zeros(trials, 1);


    % Frequency axis for reconstruction

    x = discretize(-0.5, 0.5, M);


    opt_fi.debug = 2;
    opt_fi.allow_eps_too_small = 0;
    opt_fi.nthreads = 12;

    finufft_tol = 1e-7;


    if random_signals

        % signals = randn(1, N); % + 1i*randn(trials, N);
        % signals = signals ./ vecnorm(signals, 2, 2);
        % phases = mod((0:N-1) + randn(trials, N)*s_nu, N);

        signal = randn(N, 1);
        signal = signal ./ norm(signal);

        phase = ...
            mod((0:N-1).' + randn(N, 1) * s_nu, N);

    end


    % Run multiple trials

    for i = 1:trials

        fprintf('Trial number %d\n', i);


        % Generate signal based on test mode

        if random_signals

            % signal = signals(i, :).';
            % phase = phases(i, :).';

        else

            phase = ...
                mod((0:N-1).' + randn(N, 1) * s_nu, N);

            R = randi(1) + 1;

            refl = rand(R, 1) - 0.5;
            ampl = abs(randn(R, 1));

            signal = ...
                gen_signal_from_refl( ...
                    phase, refl, ampl, true);

        end


        if compute_true

            % True reconstruction (ground truth)

            y_true = reconstruction( ...
                phase, signal, x);

        end


        tic

        y_0 = gridding_al( ...
            signal, phase, M, W, sp0, h0, alpha);

        time_0(i) = toc;


        % Gridding reconstruction using optimized kernel

        tic

        y = gridding_al( ...
            signal, phase, M, W, sp, h, alpha);

        time(i) = toc;


        tic

        y_test = gridding_al( ...
            signal, phase, M, W_test, ...
            sp_test, h_test, alpha);

        time_test(i) = toc;


        phase_fi = 2 * pi * (phase / M);

        tic

        y_fi = finufft1d1( ...
            phase_fi, signal, -1, finufft_tol, M, opt_fi);

        time_fi(i) = toc;


        if compute_true

            % Compute absolute and relative errors

            err(i, :) = abs(y_true - y);
            err_test(i, :) = abs(y_true - y_test);
            err_fi(i, :) = abs(y_true - y_fi);
            err_0(i, :) = abs(y_true - y_0);

        end

    end


    m_t = mean(time);
    m_t_test = mean(time_test);
    m_t_fi = mean(time_fi);
    m_t_0 = mean(time_0);


    if compute_true

        % Compute mean errors across trials

        m_err = mean(err, 1);
        m_err_test = mean(err_test, 1);
        m_err_fi = mean(err_fi, 1);
        m_err_0 = mean(err_0, 1);


        % Plot Mean Absolute Error (MAE)

        m = min( ...
            [min(m_err), min(m_err_fi), min(m_err_test)], ...
            [], ...
            'all') * 0.5;

        Mm = max( ...
            [max(m_err), max(m_err_fi), max(m_err_test)], ...
            [], ...
            'all') * 2;


        f_mae = figure('Name', 'MAE');

        semilogy( ...
            x, m_err_fi, ...
            'Color', green(1), ...
            'DisplayName', 'FINUFFT', ...
            'LineStyle', '-', ...
            'LineWidth', 4);

        hold on;

        semilogy( ...
            x, m_err_test, ...
            'Color', red(1), ...
            'DisplayName', 'KB ($W+1$)', ...
            'LineStyle', '-', ...
            'LineWidth', 4);

        semilogy( ...
            x, m_err, ...
            'Color', blue(0), ...
            'DisplayName', 'VO', ...
            'LineStyle', '-', ...
            'LineWidth', 4);

        semilogy( ...
            x, m_err_0, ...
            'Color', red(0), ...
            'DisplayName', 'KB ($W$)', ...
            'LineStyle', '-', ...
            'LineWidth', 4);

        % semilogy(x, sqrt(Lambda_final), ...
        %     'Color', 'black', ...
        %     'LineStyle', ':', ...
        %     'DisplayName', '$\sqrt{\Lambda C}$', ...
        %     'LineWidth', 4)

        hold off;

        title("");

        legend('Location', 'best');

        xlabel("$\omega$");
        ylabel("MAE$(\omega)$");

        fontsize(20, "points");

        ylim([m, Mm]);

    end


    ctime_vo(d) = m_t;
    ctime_kb(d) = m_t_test;
    ctime_fi(d) = m_t_fi;

end


%% ---- Timing ------------------------------------------------------------

ind = [1, 2, 3, 4, 5, 6];

fig_times = figure('Name', 'Times');

loglog( ...
    Ns(ind), ctime_fi(ind), ...
    'Color', green(1), ...
    'DisplayName', 'FINUFFT', ...
    'LineStyle', ':', ...
    'LineWidth', 2.5, ...
    'Marker', '+', ...
    'MarkerSize', 15);

hold on;

loglog( ...
    Ns(ind), ctime_kb(ind), ...
    'Color', red(1), ...
    'DisplayName', 'KB ($W+1$)', ...
    'LineStyle', ':', ...
    'LineWidth', 2.5, ...
    'Marker', '*', ...
    'MarkerSize', 15);

loglog( ...
    Ns(ind), ctime_vo(ind), ...
    'Color', blue(0), ...
    'DisplayName', 'VO', ...
    'LineStyle', ':', ...
    'LineWidth', 2.5, ...
    'Marker', 'o', ...
    'MarkerSize', 15);

hold off;

xlabel("$N$");
ylabel("time (s)");

legend('Location', 'best');

fontsize(30, "points");

ylim([1e-4, 2e0]);