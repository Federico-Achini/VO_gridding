function y = reconstruction(phase, signal, x)
%RECONSTRUCTION Compute the non-uniform Fourier reconstruction of a signal.
%
%   y = reconstruction(phase, signal, x)
%
%   Inputs
%   ------
%   phase : vector
%       Vector of phase locations (frequencies), assumed to be
%       non-uniformly spaced.
%
%   signal : vector
%       Signal values corresponding to the phase locations.
%
%   x : vector
%       Uniform grid points where the reconstruction is evaluated.
%
%   Output
%   ------
%   y : vector
%       Reconstructed signal evaluated at points x.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Input orientation

    % Ensure input vectors have correct orientations for matrix multiplication.
    phase  = orientation(phase,  'h');   % Row vector
    signal = orientation(signal, 'v');   % Column vector
    x      = orientation(x,      'v');   % Column vector


    %% Complex exponential kernel

    % Each element K(m,n) = exp(-2i * pi * x(m) * phase(n)).
    K = exp(-2i * pi * x .* phase);


    %% Reconstruction

    % Compute reconstruction by summing over signal weighted by kernel.
    y = K * signal;

end