function [error] = funRegister(p, image1, image2, scalingParam)
% p is the vector containing the deformation parameters (x translation, y
% translation).
tx = scalingParam * p(1);
ty = scalingParam * p(2);

% Compute new grid (xi,yi) using p.
n = length(image2);
[X, Y] = meshgrid(1:n, 1:n);
Xi = X + tx;
Yi = Y + ty;

% Interpolate new image.
translatedImage2 = interp2(X, Y, image2, Xi, Yi, 'linear', 0);

% Compute similarity measure between images and return as error. We take
% the negative of the error, because we are trying to minimize using
% fminsearch() but we want to maximum the normalized correlation
% coefficient.
error = -myNCC(image1, translatedImage2);
end