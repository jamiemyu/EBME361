function R = myNCC(image1, image2)
% Compute the normalized cross correlation between two input images.
% Parameters image1 and image2 are the two input images.
% Output R is the measure of how much similarity there is between two
% images (R = 1 is perfectly similar; R = 0 is perfectly dissimilar).

% Find the size of the matrices.
[row, col] = size(image1);

% Calculate the mean of each image matrix.
meanImage1 = mean2(image1);
meanImage2 = mean2(image2);

% Initialize sum variables.
num = 0;
stdImage1 = 0;
stdImage2 = 0;
for i = 1:row
    for j = 1:col
        % Compute the numerator.
        num = num + ((image1(i,j) - meanImage1) .* ...
                     (image2(i,j) - meanImage2));
        % Compute the standard deviations for each image.         
        stdImage1 = stdImage1 + (image1(i,j) - meanImage1)^2;
        stdImage2 = stdImage2 + (image2(i,j) - meanImage2)^2;
    end
end

% Compute the denominator.
denom = sqrt(stdImage1) * sqrt(stdImage2);

R = num / denom;
end