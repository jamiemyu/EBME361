function [H, filteredImageDFT, filteredImage] = gaussLPF(image, sigma)
    % Get image dimensions    
    [M, ~] = size(image);
    
    % Transform image to Fourier Domain.
    imageDFT = fft2(image);
    
    % Construct the gaussian filter of same size and sigma as f_c. 
    H = fspecial('gaussian', M, sigma);
    H = mat2gray(H);
    
    % Apply the filter by multiplying with the image in Fourier Domain
    % (which is equal to convolving in the Time Domain).
    filteredImageDFT = fftshift(H) .* imageDFT;
    
    % Transform the image back into the Time Domain.
    filteredImage = real(ifft2(filteredImageDFT));
end