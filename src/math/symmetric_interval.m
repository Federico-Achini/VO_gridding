function y = symmetric_interval(a, b, dt)
%SYMMETRIC_INTERVAL Generate a symmetric vector over [a, b].
%
%   y = symmetric_interval(a, b, dt)
%
%   Inputs
%   ------
%   a : scalar
%       Left endpoint of the interval.
%
%   b : scalar
%       Right endpoint of the interval.
%
%   dt : scalar
%       Step size.
%
%   Output
%   ------
%   y : row vector
%       Values in [a, b], symmetric around the center, excluding a
%       repeated center or endpoint when needed.
%
%   Example
%   -------
%   symmetric_interval(-2, 2, 1) -> [-2, -1, 0, 1]
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Construct symmetric interval

    half = (b + a) / 2;                 % Compute midpoint of interval
    right = half : dt : b;               % Right half starting from center
    left = flip(half - dt : -dt : a);    % Left half, symmetric to right

    if right(end) == b                   % Avoid duplicating the endpoint
        right(end) = [];
    end

    y = [left, right];                   % Combine left and right halves

end