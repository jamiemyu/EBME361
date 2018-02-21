%% Problem 1. 
clear; close all; clc
% Apply each of the noise smoothing filters (3x3 averaging, 9x9 averaging,
% 3x3 median, wiener2) to the test images, Edges_gnoise (edges + 
% Gaussian noise) and Edges_spnoise (edges + salt & pepper noise).

% Load the original noisy images.
gaussEdges = imread('Edges_gnoise.tif');
spEdges = imread('Edges_spnoise.tif');
fig1 = figure(1);
subplot(1,2,1);
imshow(gaussEdges); title('Original Gaussian Noisy Image');
subplot(1,2,2); 
imshow(spEdges); title('Original Salt & Pepper Noisy Image');
saveas(fig1, 'problem1_fig1.jpg');

% Create the averaging filters.
h1 = fspecial('average', 3); % 3x3
h2 = fspecial('average', 9); % 9x9

% Filter each image using both average filters.
gaussOut1 = imfilter(gaussEdges(:,:,1), h1, 'conv'); % 3x3
gaussOut2 = imfilter(gaussEdges(:,:,1), h2, 'conv'); % 9x9

spOut1 = imfilter(spEdges(:,:,1), h1, 'conv'); % 3x3
spOut2 = imfilter(spEdges(:,:,1), h2, 'conv'); % 9x9

% Filter each image using a median filter.
gaussOut3 = medfilt2(gaussEdges(:,:,1), [3 3]);
spOut3 = medfilt2(spEdges(:,:,1), [3 3]);

% Create the wiener2 filter. Default neighborhood of size 3x3 is chosen.
gaussOut4 = wiener2(gaussEdges(:,:,1)); 
spOut4 = wiener2(spEdges(:,:,1));

% Display the results of filtering on each noisy image.
fig2 = figure(2);
subplot(2,2,1);
imshow(gaussOut1); title('Gaussian w/ 3x3 Averaging Filter');
subplot(2,2,2);
imshow(gaussOut2); title('Gaussian w/ 9x9 Averaging Filter');
subplot(2,2,3);
imshow(gaussOut3); title('Gaussian w/ 3x3 Median Filter');
subplot(2,2,4);
imshow(gaussOut4); title('Gaussian w/ wiener2() filter');
saveas(fig2, 'problem1_fig2.jpg');

fig3 = figure(3);
subplot(2,2,1);
imshow(spOut1); title('Salt & Pepper w/ 3x3 Averaging Filter');
subplot(2,2,2);
imshow(spOut2); title('Salt & Pepper w/ 9x9 Averaging Filter');
subplot(2,2,3);
imshow(spOut3); title('Salt & Pepper w/ 3x3 Median Filter');
subplot(2,2,4);
imshow(spOut4); title('Salt & Pepper w/ wiener2() filter');
saveas(fig3, 'problem1_fig3.jpg');

%% Problem 2.
close all; clear; clc;
 
% Create a 2D image containing only Gaussian noise with a standard
% deviation of 20 gray-scale values.
gaussNoise = 20 * randn(400);
gaussImage = 150 * ones(400) + gaussNoise;
 
fig4 = figure(4);
imshow(uint8(gaussImage)); title('2D Image of Gaussian Noise w/ STD = 20');
saveas(fig4, 'problem2_fig1.jpg');
 
% Compute the standard deviation of the noise before filtering.
stdInput = std2(gaussImage(:));
 
% Create a 5x5 averaging filter.
h1 = fspecial('average', 5);
 
% Filter the image using the averaging filter.
out1 = imfilter(gaussImage, h1, 'conv');
 
% Display the image.
fig5 = figure(5);
imshow(uint8(out1)); title('Image filtered with 5x5 averaging filter');
saveas(fig5, 'problem2_fig2.jpg');
 
% Compute the standard deviation of the noise after filtering.
stdAveragingOutput = std2(out1(:));
 
% Filter the image using a median filter.
out2 = medfilt2(gaussImage, [5 5]);
 
% Display the image. 
fig6 = figure(6);
imshow(uint8(out2)); title('Image filtered with 5x5 median filter');
saveas(fig6, 'problem2_fig3.jpg');
 
% Compute the standard deviation of the noise after filtering.
stdMedianOutput = std2(out2(:));

%% Problem 3.
close all; clear; clc;

% Load the original noisy image.
spEdges = imread('Edges_spnoise.tif');
fig7 = figure(7);
imshow(spEdges); title('Original Salt & Pepper Noisy Image');
saveas(fig7, 'problem3_fig1.jpg');

% Create a 5x5 median filter (cross-shaped).
xMedFilter = [0 0 1 0 0; ...
              0 0 1 0 0; ...
              1 1 1 1 1; ...
              0 0 1 0 0; ...
              0 0 1 0 0];

% Filter the image.
medOut1 = ordfilt2(spEdges(:,:,1), ceil(nnz(xMedFilter)/2), xMedFilter);

% Create a 3x3 median filter.
medOut2 = medfilt2(spEdges(:,:,1), [3 3]);

fig8 = figure(8);
subplot(1,2,1); imshow(medOut1); title('5x5 Cross-Shaped Median Filter');
subplot(1,2,2); imshow(medOut2); title('3x3 Median Filter');
saveas(fig8, 'problem3_fig2.jpg');

stdOrigImage = std2(spEdges(:));
stdCrossImage = std(medOut1(:));
stdMedImage = std(medOut2(:));

fprintf('The standard deviation of the original image is %.4f\n', stdOrigImage);
fprintf('The standard deviation of the 5x5 cross-shape median filtered image is %.4f\n', stdCrossImage);
fprintf('The standard deviation of the 3x3 median filtered image is %.4f\n', stdMedImage);

%% Problem 4.
close all; clear; clc;

% Create a 2D image containing only Gaussian noise with a standard
% deviation of 20 gray-scale values.
gaussNoise = 20 * randn(400);
gaussImage = 150 * ones(400) + gaussNoise;
fig8 = figure(8);
imshow(uint8(gaussImage)); title('2D Image of Gaussian Noise w/ STD = 20');
saveas(fig8, 'problem4_fig1.jpg');

% Filter with kernel size 5x5 and sigma1 = 0.5, sigma2 = 5.
% Estimate three standard deviations for the filter size.
size1 = (1+ceil(0.5 * 3));
size2 = ceil(5 * 3);
gaussFilt1 = fspecial('gaussian', [size1 size1], 0.5);
gaussFilt2 = fspecial('gaussian', [size2 size2], 5);
integralFilt1 = trapz(gaussFilt1(:));
integralFilt2 = trapz(gaussFilt2(:));
fprintf('Filter with sigma = 0.5 integrates to %.4f\n', integralFilt1);
fprintf('Filter with sigma = 5 integrates to %.4f\n\n', integralFilt2);

gaussOut1 = imfilter(gaussImage, gaussFilt1, 'conv');
gaussOut2 = imfilter(gaussImage, gaussFilt2, 'conv');

% Create an averaging filter of the same kernel order.
avgFilt1 = 1 / (size1 ^ 2) * ones(size1);
avgFilt2 = 1 / (size2 ^ 2) * ones(size2);

avgOut1 = imfilter(gaussImage, avgFilt1, 'conv');
avgOut2 = imfilter(gaussImage, avgFilt2, 'conv');

% Show the results.
fig9 = figure(9);
subplot(2,2,1); imshow(uint8(gaussOut1)); title('Filtered image, \sigma = 0.5');
subplot(2,2,2); imshow(uint8(gaussOut2)); title('Filtered image, \sigma = 5');
subplot(2,2,3); imshow(uint8(avgOut1)); title('Average filtered image for order of \sigma = 0.5');
subplot(2,2,4); imshow(uint8(avgOut2)); title('Average filtered image for order of \sigma = 5');
saveas(fig9, 'problem4_fig2.jpg');

% Measure the DC gain of the Gaussian kernels.
dcGain1 = sum(sum(gaussFilt1));
dcGain2 = sum(sum(gaussFilt2));
fprintf('DC gain (for sigma = 0.5) = %.4f\n', dcGain1);
fprintf('DC gain (for sigma = 5) = %.4f\n\n', dcGain2);

% Measure the mean of input/output noise.
meanIn = mean(gaussImage(:));
meanOut1 = mean(gaussOut1(:));
meanOut2 = mean(gaussOut2(:));

% Measure the standard deviation of input/output noise.
stdIn = std2(gaussImage(:));
stdOut1 = std2(gaussOut1(:));
stdOut2 = std2(gaussOut2(:));

% Normalize kernels.
normInput = (gaussImage - min(gaussImage(:))) / ...
            (max(gaussImage(:)) - min(gaussImage(:)));
normGaussOut1 = (gaussOut1 - min(gaussOut1(:))) / ...
    (max(gaussOut1(:)) - min(gaussOut1(:)));
normGaussOut2 = (gaussOut2 - min(gaussOut2(:))) / ...
    (max(gaussOut2(:)) - min(gaussOut2(:)));

% Measure the mean of normalized input/output noise.
normMeanIn = mean(normInput(:));
normMeanOut1 = mean(normGaussOut1(:));
normMeanOut2 = mean(normGaussOut2(:));

% Measure the std of normalized input/output noise.
normStdIn = std2(normInput(:));
normStdOut1 = std2(normGaussOut1(:));
normStdOut2 = std2(normGaussOut2(:));

% Measure the average filtered mean noise.
meanInAvg1 = mean(avgOut1(:));
meanInAvg2 = mean(avgOut2(:));

% Measure the average filtered standard deviation of noise.
stdInAvg1 = std2(avgOut1(:));
stdInAvg2 = std2(avgOut2(:));

fprintf('Input mean = %.4f, Output 1 mean = %.4f, Output 2 mean = %.4f\n', meanIn, meanOut1, meanOut2);
fprintf('Norm input mean = %.4f, Norm output 1 mean = %.4f, Norm output 2 mean = %.4f\n\n', normMeanIn, normMeanOut1, normMeanOut2);
fprintf('Input std = %.4f, Output 1 std = %.4f, Output 2 std = %.4f\n', stdIn, stdOut1, stdOut2);
fprintf('Norm input std = %.4f, Norm output 1 std = %.4f, Norm output 2 std = %.4f\n\n', normStdIn, normStdOut1, normStdOut2);
fprintf('Average Filtered Output 1 mean = %.4f, Average Filtered Output 2 mean = %.4f\n', meanInAvg1, meanInAvg2);
fprintf('Average Filtered Output 1 std = %.4f, Average Filtered Output 2 std = %.4f\n', stdInAvg1, stdInAvg2);
