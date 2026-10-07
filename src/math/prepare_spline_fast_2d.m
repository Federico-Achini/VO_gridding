function sk = prepare_spline_fast_2d(spl1, spl2, W, D)
%PREPARE_SPLINE_FAST_2D Prepare spline data for 2-D kernel evaluation.
%
%   sk = prepare_spline_fast_2d(spl1, spl2, W, D)
%
%   Inputs
%   ------
%   spl1, spl2 : spline pairs
%       Spline pairs (real, imaginary) for the kernel factor along
%       dimension 1 and dimension 2.
%
%   W : scalar or 1x2 vector
%       Kernel half-width. A scalar value is used for both dimensions.
%
%   D : scalar or 1x2 vector
%       Discretization factor. A scalar value is used for both dimensions.
%
%   Output
%   ------
%   sk : 1x2 struct array
%       Prepared spline data for the two dimensions.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Expand scalar parameters

    if isscalar(W)
        W = [W W];
    end

    if isscalar(D)
        D = [D D];
    end


    %% Validate parameters

    if numel(W) ~= 2 || numel(D) ~= 2

        error('W and D must each be a scalar or a 1x2 vector.');

    end


    %% Prepare spline data

    sk = [prepare_spline_fast(spl1, W(1), D(1)), ...
          prepare_spline_fast(spl2, W(2), D(2))];   % 1x2 struct array

end