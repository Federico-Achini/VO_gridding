function y = crop_norm(img, x, y, xcenter, ycenter, xwidth, ywidth)
%CROP_NORM Compute the norm of a cropped region of an image.
%
%   y = crop_norm(img, x, y, xcenter, ycenter, xwidth, ywidth)
%
%   Inputs
%   ------
%   img : array
%       Input image.
%
%   x, y : vectors
%       Spatial coordinates corresponding to the image dimensions.
%
%   xcenter, ycenter : scalar
%       Center coordinates of the cropped region.
%
%   xwidth, ywidth : scalar
%       Widths of the cropped region along the x and y dimensions.
%
%   Output
%   ------
%   y : scalar
%       Norm of the cropped image region.
%
%   Author: Federico Achini
%   Date:   August 2025
%
%   ----------------------------------------------------------------------

    %% Crop indices

    ix = abs(x - xcenter) < xwidth/2;
    iy = abs(y - ycenter) < ywidth/2;


    %% Compute norm

    y = norm(img(ix, iy));

end