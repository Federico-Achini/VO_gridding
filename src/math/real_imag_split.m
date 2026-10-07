function c = real_imag_split(x)
%REAL_IMAG_SPLIT Split a complex vector into real and imaginary parts.
%
%   c = real_imag_split(x)
%
%   Input
%   -----
%   x : array
%       Complex-valued vector or array.
%
%   Output
%   ------
%   c : matrix
%       Real-valued matrix with two columns:
%           - First column contains the real parts of x.
%           - Second column contains the imaginary parts of x.
%
%   Description
%   -----------
%   This function separates complex data into real and imaginary
%   components for further processing or storage.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    c = [real(x), imag(x)];

end