function four = gridding_al(signal, phase, M, W, sk, h, alpha)
%GRIDDING_AL Perform 1-D gridding with a precomputed kernel.
%
%   four = gridding_al(signal, phase, M, W, sk, h, alpha)
%
%   Inputs
%   ------
%   signal : column vector
%       Input signal samples.
%
%   phase : column vector
%       Phase coordinates.
%
%   M : scalar
%       Output grid size.
%
%   W : scalar
%       Kernel width.
%
%   sk : struct
%       Kernel data.
%
%   h : column vector
%       Correction function.
%
%   alpha : scalar
%       Oversampling factor.
%
%   Output
%   ------
%   four : column vector
%       Gridded and corrected Fourier-domain data.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        signal (:, 1) double
        phase  (:, 1) double
        M      (1, 1) double
        W      (1, 1) double
        sk     (1, 1) struct
        h      (:, 1) double
        alpha  (1, 1) double
    end


    %% Input orientation

    signal = signal(:);
    phase  = phase(:);
    h      = h(:);


    %% Oversampled grid

    aM = floor(M * alpha);


    %% Kernel gridding

    % Kernel eval + alpha*phase + grid indexing + scatter-accumulate,
    % all fused in one MEX call

    %X = grid_accumulate_mex(sk, phase, alpha, signal, aM, W);

    X = grid_par(sk, phase, alpha, signal, aM, W);


    %% FFT

    four = fftshift(fft(X));


    %% Crop

    inds = crop_oversampling(M, alpha);
    four = four(inds);


    %% Correction

    four = h .* four;

end