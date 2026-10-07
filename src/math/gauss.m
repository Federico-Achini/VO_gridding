function y = gauss(x, mu, sigma)
%GAUSS Evaluate a Gaussian function.
%
%   y = gauss(x, mu, sigma)
%
%   Inputs
%   ------
%   x : array
%       Evaluation points.
%
%   mu : scalar
%       Center of the Gaussian.
%
%   sigma : scalar
%       Width parameter of the Gaussian.
%
%   Output
%   ------
%   y : array
%       Gaussian function evaluated at x.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    y = exp(-(x - mu).^2 / sigma);

end