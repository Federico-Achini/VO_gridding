function four = gridding_al_2d(signal, phase, M, W, sk, h, alpha)
%GRIDDING_AL_2D Perform 2-D gridding with a separable kernel.
%
%   four = gridding_al_2d(signal, phase, M, W, sk, h, alpha)
%
%   Inputs
%   ------
%   signal : N x 1 double
%       Input signal samples.
%
%   phase : N x 2 double
%       Phase coordinates. Column 1 corresponds to dimension 1;
%       column 2 corresponds to dimension 2.
%
%   M : scalar or 1 x 2 double
%       Output grid size.
%
%   W : scalar or 1 x 2 double
%       Kernel width.
%
%   sk : 1 x 2 struct array
%       Kernel data from prepare_spline_fast_2d.
%
%   h : array or cell array
%       Correction matrix of size M(1) x M(2), or {h1, h2} of
%       column vectors. For a separable kernel, the correction is
%       the outer product h1*h2.'.
%
%   alpha : scalar or 1 x 2 double
%       Oversampling factor.
%
%   Output
%   ------
%   four : matrix
%       Gridded and corrected Fourier-domain data.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        signal (:, 1) double
        phase  (:, 2) double
        M      (1, :) double
        W      (1, :) double
        sk     (1, :) struct
        h
        alpha  (1, :) double
    end


    %% Expand scalar parameters

    M     = expand2(M,     'M');
    W     = expand2(W,     'W');
    alpha = expand2(alpha, 'alpha');


    %% Correction

    if iscell(h)
        h = h{1}(:) .* h{2}(:).';
    end


    %% Oversampled grid

    aM = floor(M .* alpha);                          % 1x2 oversampled grid size

    % Kernel eval + alpha*phase + grid indexing + scatter-accumulate, fused
    X = grid_par2d(sk, phase, alpha, signal, aM, W); % aM(1) x aM(2)


    %% FFT

    four = fftshift(fft2(X));


    %% Crop

    four = four( ...
        crop_oversampling(M(1), alpha(1)), ...
        crop_oversampling(M(2), alpha(2)));


    %% Correction

    if ~isequal(size(h), size(four))
        error( ...
            'gridding_al_2d:hsize', ...
            'h must be %d-by-%d.', ...
            size(four, 1), ...
            size(four, 2));
    end

    four = h .* four;

end


function v = expand2(v, name)

    if isscalar(v)
        v = [v v];

    elseif numel(v) ~= 2
        error( ...
            'gridding_al_2d:size', ...
            '%s must be a scalar or a 1x2 vector.', ...
            name);
    end

end