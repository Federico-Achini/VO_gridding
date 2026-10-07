function y = orientation(x, how)
%ORIENTATION Reshape a vector to a specified orientation.
%
%   y = orientation(x, how)
%
%   Inputs
%   ------
%   x : vector
%       Input vector, either row or column oriented.
%
%   how : character vector
%       Desired orientation:
%           'v' - column vector (N x 1)
%           'h' - row vector (1 x N)
%
%   Output
%   ------
%   y : vector
%       Reshaped vector with the specified orientation.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    switch how

        case 'v'  % Convert to column vector
            y = reshape(x, [], 1);

        case 'h'  % Convert to row vector
            y = reshape(x, 1, []);

        otherwise
            error('Orientation type "%s" not available.', how);

    end

end