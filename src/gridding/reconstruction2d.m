function img = reconstruction2d(phase, signal, Mx, My)
%RECONSTRUCTION2D Reconstruct a 2-D image from nonuniform frequency samples.
%
%   img = reconstruction2d(phase, signal, Mx, My)
%
%   Inputs
%   ------
%   phase : array
%       Nonuniform frequency coordinates.
%
%   signal : array
%       Signal samples corresponding to the frequency coordinates.
%
%   Mx : scalar
%       Number of reconstruction points along the x dimension.
%
%   My : scalar
%       Number of reconstruction points along the y dimension.
%
%   Output
%   ------
%   img : My x Mx array
%       Reconstructed 2-D image.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Reconstruction grid

    x = discretize(-0.5, 0.5, Mx);   % [1 x Mx]
    y = discretize(-0.5, 0.5, My);   % [1 x My]


    %% Frequency coordinates

    [~, Ny, Nx] = size(phase);

    % Flatten frequency coordinates to [1 x Nx*Ny].
    kx = reshape(phase(2, :, :), 1, []);   % [1 x Nx*Ny]
    ky = reshape(phase(1, :, :), 1, []);   % [1 x Nx*Ny]
    s  = reshape(signal,         1, []);   % [1 x Nx*Ny]


    %% Exponential matrices

    % Build exponential matrices — no loop, pure broadcasting.
    %
    % Ex(j, m) = exp(2pi*i * x(j) * kx(m))  -> [Mx x Nx*Ny]
    % Ey(i, m) = exp(2pi*i * y(i) * ky(m))  -> [My x Nx*Ny]

    Ex = exp(2i * pi * x(:) * kx);    % [Mx x Nx*Ny]
    Ey = exp(2i * pi * y(:) * ky);    % [My x Nx*Ny]


    %% Signal weighting

    % Weight the signal into Ex.
    % Each column of Ex is scaled by the corresponding s value.
    Ews = Ex .* s;                     % [Mx x Nx*Ny]


    %% Reconstruction

    % Contract over the frequency dimension.
    img = (Ey * Ews.') / (Mx * My);   % [My x Mx]

end