function y = eval_spline(sp, t)
%EVAL_SPLINE Fast direct cubic spline evaluation.
%
%   y = eval_spline(sp, t)
%
%   Inputs
%   ------
%   sp : spline structure
%       Spline representation and associated metadata.
%
%   t : array
%       Points at which the spline is evaluated.
%
%   Output
%   ------
%   y : array
%       Spline values evaluated at t, with the same size as t.
%
%   Description
%   -----------
%   Direct cubic spline evaluation without PPVAL or DISCRETIZE.
%   Supports matrix-valued evaluation points.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    sz = size(t);

    t = t(:);


    %% Find spline interval

    if sp.uniform

        % Direct interval calculation.
        idx = floor((t - sp.x0) * sp.invDx) + 1;

        % Handle right endpoint.
        idx(idx > sp.nPieces) = sp.nPieces;

    else

        % General non-uniform spline.
        idx = zeros(size(t));

        for k = 1:sp.nPieces

            if k == 1

                mask = t <= sp.br(2);

            else

                mask = t > sp.br(k) & t <= sp.br(k+1);

            end

            idx(mask) = k;

        end

        idx(idx == 0) = 1;

    end


    %% Local coordinate

    x = t - sp.br(idx);


    %% Direct cubic evaluation

    c = sp.cr;

    a  = c(idx,1);
    b  = c(idx,2);
    cc = c(idx,3);
    d  = c(idx,4);

    y = ((a .* x + b) .* x + cc) .* x + d;


    %% Restore shape

    y = reshape(y, sz);

end