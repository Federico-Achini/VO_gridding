function spl = real_imag_spline(t, fun)
%REAL_IMAG_SPLINE Create spline representations of real and imaginary parts.
%
%   spl = real_imag_spline(t, fun)
%
%   Inputs
%   ------
%   t : vector
%       Sample points (knots) where the function is evaluated.
%
%   fun : array
%       Complex-valued function samples at points t.
%
%   Output
%   ------
%   spl : array
%       Array containing spline structures:
%           spl(1) - spline of the real part of fun
%           spl(2) - spline of the imaginary part of fun
%
%   Description
%   -----------
%   Constructs separate spline interpolants for the real and imaginary
%   parts of a complex-valued function, enabling efficient evaluation.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Real and imaginary spline representations

    % Create spline for real part.
    spl_re = spline(t, real(fun));

    % Create spline for imaginary part.
    spl_im = spline(t, imag(fun));


    %% Combine spline representations

    % Combine into output (as an array of splines).
    spl = [spl_re, spl_im];

end