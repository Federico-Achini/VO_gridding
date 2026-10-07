function c = spline_to_legendre(spl, L, W)
%SPLINE_TO_LEGENDRE Convert a spline representation to Legendre coefficients.
%
%   c = spline_to_legendre(spl, L, W)
%
%   Inputs
%   ------
%   spl : spline representation
%       Spline representations of the real and imaginary components.
%
%   L : positive integer
%       Number of Legendre coefficients.
%
%   W : positive scalar
%       Half-width of the support interval.
%
%   Output
%   ------
%   c : L x 2 array
%       Legendre coefficients. The first column contains the real
%       coefficients and the second column contains the imaginary
%       coefficients.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    arguments
        spl
        L (1,1) double {mustBePositive, mustBeInteger}
        W (1,1) double {mustBePositive}
    end


    %% Initialize coefficients

    c = zeros(L, 2);


    %% Convert spline components

    for j = 1:2

        [br, coefs] = unmkpp(spl(j));
        D = numel(br) - 1;

        [q, w] = gauss_legendre(ceil((L + 3) / 2));

        n = (0:L-1).';
        scale = (2*n + 1) / (2*W);

        cj = zeros(L, 1, 'like', coefs);
        comp = zeros(L, 1, 'like', coefs);


        %% Integrate over spline intervals

        for d = 1:D

            xl = br(d);
            xr = br(d+1);

            xc = (xl + xr) / 2;
            h  = (xr - xl) / 2;

            x = xc + h*q;
            u = x / W;

            z = x - xl;
            f = ((coefs(d,1).*z + coefs(d,2)).*z + coefs(d,3)).*z ...
                + coefs(d,4);

            P = legendre_poly(L, u);

            contribution = h * scale .* (P * (w .* f));


            %% Kahan summation

            y = contribution - comp;
            t = cj + y;
            comp = (t - cj) - y;
            cj = t;

        end

        c(:,j) = cj;

    end

end