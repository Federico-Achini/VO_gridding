function ind = crop_oversampling(N, alpha)
%CROP_OVERSAMPLING Compute indices for cropping an oversampled FFT array.
%
%   ind = crop_oversampling(N, alpha)
%
%   Inputs
%   ------
%   N : positive integer
%       Size of the non-oversampled array.
%
%   alpha : positive scalar
%       Oversampling factor. Default: 1.
%
%   Output
%   ------
%   ind : vector
%       Indices used to crop the oversampled fftshifted array.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        N (1,1) {mustBeInteger, mustBePositive}
        alpha (1,1) {mustBePositive} = 1
    end


    %% Oversampled size

    aN = floor(alpha * N);


    %% Cropping indices

    ind = floor(aN/2) + 1 - floor(N/2) + (0:N-1);

end