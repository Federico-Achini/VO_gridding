function [c_min, hst] = optimization(fun, c, c00)
%OPTIMIZATION Perform optimization of kernel coefficients.
%
%   [c_min, hst] = optimization(fun, c, c00)
%
%   Inputs
%   ------
%   fun : function handle
%       Objective function to minimize.
%
%   c : array
%       Initial guess for the optimization variables.
%
%   c00 : array
%       Reference array used to set the TypicalX option.
%
%   Outputs
%   -------
%   c_min : array
%       Optimized coefficients.
%
%   hst : array
%       History of objective function values recorded at each iteration.
%
%   Description
%   -----------
%   This function uses MATLAB's FMINUNC with the quasi-Newton algorithm
%   to minimize the provided objective function. The optimization
%   progress is tracked through an output function, recording the
%   objective value at every iteration.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Nonlinear constraint function
    %
    % % Nonlinear constraint function:
    % % - Inequality constraints: none (empty vector)
    % % - Equality constraint: norm(x)^2 - 1 = 0 (enforces unit norm)
    % function [c, ceq, grad_c, grad_ceq] = constr(x)
    %     N = numel(x);
    %     x = resize(x, [N, 1]);
    %     c = []; % No inequality constraints
    %     ceq = [(norm(x)^2 - 1)]; % Equality constraint: norm(x) = 1
    %     grad_c = [];
    %     %grad_ceq = [];
    %     grad_ceq = [2*x];
    % end


    %% Optimization history

    % Initialize history array to store objective values.
    hst = [];


    %% Output function

    % Output function to record optimization progress.
    function stop = my_iter(x, optimValues, state)
        stop = false; % Continue optimization

        if isequal(state, 'iter')
            hst = [hst; optimValues.fval]; % Append current objective value
        end
    end


    %% Optimization options

    % Set optimization options for fmincon.
    options = optimoptions('fmincon', ...
        'Algorithm', 'interior-point', ...
        'OutputFcn', @my_iter, ...
        'Display', 'iter-detailed', ...
        "SubproblemAlgorithm", "cg", ...
        'MaxIterations', 1e3, ...
        'MaxFunctionEvaluations', 1e5, ..., ...
        'StepTolerance', 1e-8, ...
        "ConstraintTolerance", 1-1e-16, ...
        "OptimalityTolerance", 1e-8, ...
        "SpecifyConstraintGradient", true);


    options = optimoptions("fminunc", ...
        "Algorithm", "quasi-newton", ...
        "Display", "iter-detailed", ...
        'MaxIterations', 1e3, ...
        'MaxFunctionEvaluations', 1e6, ...
        'OptimalityTolerance', 1e-16, ...
        'FiniteDifferenceType', 'central', ...
        'FiniteDifferenceStepSize', 1e-9, ...
        'TypicalX', ones(size(c00)), ...
        'StepTolerance', 1e-8, ...
        'ObjectiveLimit', -inf);


    options = optimoptions('fminunc', ...
        'Algorithm', 'quasi-newton', ...
        'SpecifyObjectiveGradient', true, ...
        'Display', 'iter-detailed', ...
        'MaxIterations', 1000, ...
        'MaxFunctionEvaluations', 1e5, ...
        'OptimalityTolerance', 1e-8, ...
        'StepTolerance', 1e-10);


    %% Optimization

    % Run constrained optimization with initial guess c.
    %c_min = fmincon(fun, c, [], [], [], [], [], [], @constr, options);

    c_min = fminunc(fun, c, options);

end