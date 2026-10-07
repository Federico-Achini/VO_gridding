function guess = shifted_guess(W, D, x0, ker, mode)
%SHIFTED_GUESS Generate a shifted initial guess based on a DPSS vector.
%
%   guess = shifted_guess(W, D, x0)
%
%   Inputs
%   ------
%   W : scalar
%       Kernel half-bandwidth parameter.
%
%   D : scalar
%       Discretization factor (samples per unit interval).
%
%   x0 : scalar
%       Shift parameter (frequency shift).
%
%   ker : vector, optional
%       Kernel vector used as the initial basis. If empty, the first DPSS
%       vector is computed.
%
%   mode : character vector, optional
%       Endpoint convention passed to DISCRETIZE. Default: 'left'.
%
%   Output
%   ------
%   guess : complex vector
%       Shifted initial guess.
%
%   Description
%   -----------
%   This function generates an initial guess vector for kernel optimization
%   by modulating the first Discrete Prolate Spheroidal Sequence (DPSS)
%   with a complex exponential corresponding to a frequency shift x0.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        W
        D
        x0
        ker = []
        mode = 'left'
    end


    %% Number of samples

    N = 2*W*D;

    if strcmp(mode, 'both')
        N = N + 1;
    end


    %% Discretization

    % Discretize interval [-W, W) into N points.
    nu = discretize(-W, W, N, mode);


    %% DPSS kernel

    if isempty(ker)
        % Compute first DPSS vector of length N and bandwidth W.
        ker = dpss(N, W, 1);
    else
        assert(length(ker) == N);
    end

    ker = orientation(ker, 'v');


    %% Frequency shift

    % Apply frequency shift via complex exponential modulation.
    guess = ker .* exp(2i * pi * x0 * nu);

end