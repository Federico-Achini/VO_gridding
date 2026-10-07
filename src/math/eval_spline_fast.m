function y = eval_spline_fast(spl, t)
%EVAL_SPLINE_FAST Evaluate a uniformly spaced cubic spline efficiently.
%
%   y = eval_spline_fast(spl, t)
%
%   Inputs
%   ------
%   spl : spline structure
%       Spline representation containing the breakpoints and coefficients.
%
%   t : array
%       Points at which the spline is evaluated.
%
%   Output
%   ------
%   y : array
%       Spline values evaluated at t, with the same size as t.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Input and spline data

    sz = size(t);

    t = t(:);

    br = spl.breaks;

    cr = spl(1).coefs + 1i*spl(2).coefs;

    nPieces = size(cr, 1);

    dx    = br(2) - br(1);
    invDx = 1 / dx;
    x0    = br(1);


    %% Interval selection

    idx = floor((t - x0) .* invDx) + 1;


    %% Endpoint handling

    idx(t >= br(end)) = nPieces;


    %% Clamp

    idx = max(idx, 1);
    idx = min(idx, nPieces);


    %% Local coordinate

    u = t - (x0 + (idx - 1).*dx);


    %% Cubic Horner evaluation

    y = ((cr(idx,1).*u + cr(idx,2)).*u + cr(idx,3)).*u + cr(idx,4);

    y = reshape(y, sz);

end