function [H, filteredImageDFT, filteredImage] = idealLPF(image,f_c)
    % Get image dimensions    
    [M, N] = size(image);
    
    % Transform image to Fourier Domain.
    imageDFT = fft2(double(image));
    
    % Initialize new set of point coordinates.
    u = 0:(M-1);
    v = 0:(N-1);
    
    % Find indices of x- and y- coordinates outside half the dimensions.
    idx = find(u > M/2);
    idy = find(v > N/2);
    u(idx) = u(idx) - M;
    v(idy) = v(idy) - N;
    
    % Create a grid from the new coordinates.
    [V,U] = meshgrid(v,u);
    
    % Use the formula of a circle to construct the logical bounds of the
    % filter.
    D = sqrt(U.^2+V.^2);
    
    % Construct the filter by locating points that are within the radius of
    % the cutoff frequency.
    H = double(D <= f_c);
    
    % Apply the filter by multiplying with the image in Fourier Domain
    % (which is equal to convolving in the Time Domain).
    filteredImageDFT = H .* imageDFT;
    
    % Transform the image back into the Time Domain.
    filteredImage = real(ifft2(double(filteredImageDFT)));
end