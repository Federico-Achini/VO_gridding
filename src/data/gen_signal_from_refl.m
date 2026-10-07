function signal = gen_signal_from_refl(phase, refl, ampl, realBasis)
%GEN_SIGNAL_FROM_REFL Generate a normalized signal from reflectors.
%
%   signal = gen_signal_from_refl(phase, refl, ampl, realBasis)
%
%   Inputs
%   ------
%   phase : vector
%       Vector of phase values.
%
%   refl : vector
%       Vector of reflector frequencies/delays.
%
%   ampl : vector
%       Reflector amplitudes. Empty => random amplitudes.
%
%   realBasis : logical
%       true  - cosine basis.
%       false - complex exponential basis.
%
%   Output
%   ------
%   signal : column vector
%       Normalized signal with unit 2-norm.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Input orientation

    % Force orientations without calling orientation(), avoiding its overhead.
    phase = phase(:).';   % 1 x N
    refl  = refl(:);      % K x 1


    %% Amplitudes

    if isempty(ampl)
        ampl = rand(size(refl));
    else
        ampl = ampl(:);

        if numel(ampl) ~= numel(refl)
            error("Reflectivity and amplitudes must have the same number of elements.");
        end
    end


    %% Signal dimensions

    N = numel(phase);
    K = numel(refl);

    if N == 0
        signal = zeros(0, 1);
        return;
    end

    if K == 0
        signal = zeros(N, 1);
        return;
    end


    %% Block size

    % Chunk reflectors to avoid allocating a K x N matrix when K is large.
    % Increase/decrease this depending on memory and CPU characteristics.
    targetElements = 2e5;
    blockSize = max(1, floor(targetElements / max(N, 1)));


    %% Signal construction

    if realBasis
        signal = zeros(1, N);

        for k0 = 1:blockSize:K
            k1 = min(k0 + blockSize - 1, K);

            % Kblock x N temporary, substantially smaller than K x N.
            arg = 2*pi * (refl(k0:k1) .* phase);

            % Sum immediately so the full components matrix is never retained.
            signal = signal + sum(ampl(k0:k1) .* cos(arg), 1);
        end

    else
        signal = complex(zeros(1, N));

        for k0 = 1:blockSize:K
            k1 = min(k0 + blockSize - 1, K);

            arg = 2*pi * (refl(k0:k1) .* phase);

            signal = signal + sum(ampl(k0:k1) .* exp(1i * arg), 1);
        end
    end


    %% Normalize

    % Normalize without creating another unnecessary copy.
    signal = signal(:);
    signal = signal / norm(signal);

end