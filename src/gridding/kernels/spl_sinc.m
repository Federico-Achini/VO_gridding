function spl = spl_sinc(W, D, x0)
%SPL_SINC Construct a spline representation from a sampled sinc kernel.
%
%   spl = spl_sinc(W, D, x0)
%
%   Inputs
%   ------
%   W : scalar
%       Kernel width.
%
%   D : scalar
%       Number of samples per unit interval.
%
%   x0 : scalar
%       Frequency or shift parameter.
%
%   Output
%   ------
%   spl : structure
%       Spline representation of the sampled kernel.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Kernel grid

    nu = discretize(-W, W, 2*D*W);
    x = 0:2*W*D-1;


    %% Sampled signal

    u = exp(2i*pi*x0*x);


    %% Fourier transform

    k = fftshift(fft(u));


    %% Spline representation

    spl = real_imag_spline(nu, k);

end