%% Problem 1

% No coding needed here.

%% Problem 2
close all; clear; clc;

% Load the images.
wheel = imread('wheel.tif');
edges = rgb2gray(imread('edges.tif'));
[m1, n1] = size(wheel);

wheelFT = fft2(double(wheel));
edgesFT = fft2(double(edges));

% a. Apply an ideal lpf at a frequency which is 1/3 the maximum frequency.
% Maximum frequency occurs at +/- N/2 in each dir. of an unshifted NxN DFT.
magWheelFT = abs(wheelFT);
magEdgesFT = abs(edgesFT);
maxFreqWheel = magWheelFT(m1/2, n1/2);
maxFreqEdges = magEdgesFT(m1/2, n1/2);

% Compute cutoff frequencies as 1/3 the maximum frequencies.
fcWheel = 1/3 * maxFreqWheel;
fcEdges= 1/3 * maxFreqEdges;

% Construct an ideal low-pass filter which is a circle shape of 0s and 1s.
[~, lpfWheelDFT, lpfWheel] = idealLPF(wheel, fcWheel);
[edgesH, lpfEdgesDFT, lpfEdges] = idealLPF(edges, fcEdges);

% b. Apply a Gaussian lpf with a 1/2 width at half maximum amplitude which is
% 1/3 the maximum frequency.
[~, gaussWheelDFT, gaussWheel] = gaussLPF(wheel, fcWheel);
[gaussEdgesH, gaussEdgesDFT, gaussEdges] = gaussLPF(edges, fcEdges);

% c. Print out these images for the report.
% i. Original images in the cartesian domain.
fig1 = figure(1);
subplot(1,2,1); imshow(wheel, []); title('Original Wheel Image');
subplot(1,2,2); imshow(edges, []); title('Original Edges Image');
saveas(fig1, 'fig1_hw6.jpg');

% ii. Original images in the frequency domain.
fig2 = figure(2);
subplot(1,2,1); imshow(log(abs(fftshift(wheelFT))), []); title('Wheel Image in FD');
subplot(1,2,2); imshow(log(abs(fftshift(edgesFT))), []); title('Edges Image in FD');
saveas(fig2, 'fig2_hw6.jpg');

% iii. Images with Ideal LPF in cartesian domain.
fig3 = figure(3);
subplot(1,2,1); imshow(lpfWheel, []); title('Ideal LPF Wheel Image');
subplot(1,2,2); imshow(lpfEdges, []); title('Ideal LPF Edges Image');
saveas(fig3, 'fig3_hw6.jpg');

% iv. Images with Ideal LPF in frequency domain.
fig4 = figure(4);
subplot(1,2,1); imshow(log(abs(fftshift(lpfWheelDFT))), []); title('Ideal LPF Wheel Image in FD');
subplot(1,2,2); imshow(log(abs(fftshift(lpfEdgesDFT))), []); title('Ideal LPF Edges Image in FD');
saveas(fig4, 'fig4_hw6.jpg');

% v. Images with Gaussian LPF in cartesian domain.
fig5 = figure(5);
subplot(1,2,1); imshow(gaussWheel, []); title('Gauss LPF Wheel Image');
subplot(1,2,2); imshow(gaussEdges, []); title('Gauss LPF Edges Image');
saveas(fig5, 'fig5_hw6.jpg');

% vi. Images with Gaussian LPF in frequency domain.
fig6 = figure(6);
subplot(1,2,1); imshow(log(abs(fftshift(gaussWheelDFT))), []); title('Gauss LPF Wheel Image in FD');
subplot(1,2,2); imshow(log(abs(fftshift(gaussEdgesDFT))), []); title('Gauss LPF Edges Image in FD');
saveas(fig6, 'fig6_hw6.jpg');

% vii. Ideal LPF in frequency domain.
fig7 = figure(7);
imshow(fftshift(edgesH)); title('Ideal Low-Pass Filter');
saveas(fig7, 'fig7_hw6.jpg');

% viii. Gaussian LPF in frequency domain.
fig8 = figure(8);
imshow(gaussEdgesH); title('Gaussian Low-Pass Filter');
saveas(fig8, 'fig8_hw6.jpg');

% d. No code.
% e. No code.

%% f. Apply an ideal LPF and a Gaussian LPF to an image consisting of a
% single white dot on a black background.

% Create a 256x256 black background image.
dotImg = zeros(50, 50);

% Set the center pixel to a white dot.
dotImg(50/2, 50/2) = 1;

fig0 = figure(13);
subplot(1,3,1); imshow(dotImg); title('2D Impulse Function Image');

% (i) Apply an ideal LPF (select a cutoff frequency of 20)
fc = 15;
[~, ~, lpfDotImg] = idealLPF(dotImg, fc);
figure(13);
subplot(1,3,2); imshow(lpfDotImg); title('2D Impulse Function with Ideal LPF');

% (ii) Apply a Gaussian LPF.
[~, ~, gaussDotImg] = gaussLPF(dotImg, fc);
figure(13)
subplot(1,3,3); imshow(gaussDotImg); title('2D Impulse Function with Gaussian LPF');
saveas(fig0, 'fig13_hw6.jpg');

%% Problem 3
close all; clear; clc;

lena = imread('lena.bmp');
iris = imread('iris-illustration.bmp');

% a. Convert images to grayscale format.
lena = im2double(rgb2gray(lena));
iris = im2double(rgb2gray(iris));

% Calculate DFT and enlarge the spectral representation with zero padding.
lenaDFT = fft2(lena, 2 * size(lena, 1), 2 * size(lena, 2));
irisDFT = fft2(iris, 2 * size(iris, 1), 2 * size(iris, 2));

% b. Plot the log magnitude of the 2D DFT of the grayscale image, with
% center shifted.
fig9 = figure(9);
subplot(1,2,1); imshow(log(abs(fftshift(lenaDFT))), []); title('Log of Magnitude of Lena Image');
subplot(1,2,2); imshow(log(abs(fftshift(irisDFT))), []); title('Log of Magnitude of Iris Image');
saveas(fig9, 'hw6_fig9.jpg');

%% c. Apply truncation windows to keep 30%, 15%, and 5% of the DFT
% coefficients.

% Find the size of the original DFT.
[m1, n1] = size(lenaDFT);
[m2, n2] = size(irisDFT);

% Lena x- and y-coordinates
xMinLena30 = ceil(m1/2 - m1 * sqrt(0.3 / 2));
xMaxLena30 = ceil(m1/2 + m1 * sqrt(0.3 / 2));
yMinLena30 = ceil(n1/2 - n1 * sqrt(0.3 / 2));
yMaxLena30 = ceil(n1/2 + n1 * sqrt(0.3 / 2));
xMinLena15 = ceil(m1/2 - m1 * sqrt(0.15 / 2));
xMaxLena15 = ceil(m1/2 + m1 * sqrt(0.15 / 2));
yMinLena15 = ceil(n1/2 - n1 * sqrt(0.15 / 2));
yMaxLena15 = ceil(n1/2 + n1 * sqrt(0.15 / 2));
xMinLena5 = ceil(m1/2 - m1 * sqrt(0.05 / 2));
xMaxLena5 = ceil(m1/2 + m1 * sqrt(0.05 / 2));
yMinLena5 = ceil(n1/2 - n1 * sqrt(0.05 / 2));
yMaxLena5 = ceil(n1/2 + n1 * sqrt(0.05 / 2));

% Iris x- and y-coordinates
xMinIris30 = ceil(m2/2 - m2 * sqrt(0.3 / 2));
xMaxIris30 = ceil(m2/2 + m2 * sqrt(0.3 / 2));
yMinIris30 = ceil(n2/2 - n2 * sqrt(0.3 / 2));
yMaxIris30 = ceil(n2/2 + n2 * sqrt(0.3 / 2));
xMinIris15 = ceil(m2/2 - m2 * sqrt(0.15 / 2));
xMaxIris15 = ceil(m2/2 + m2 * sqrt(0.15 / 2));
yMinIris15 = ceil(n2/2 - n2 * sqrt(0.15 / 2));
yMaxIris15 = ceil(n2/2 + n2 * sqrt(0.15 / 2));
xMinIris5 = ceil(m2/2 - m2 * sqrt(0.05 / 2));
xMaxIris5 = ceil(m2/2 + m2 * sqrt(0.05 / 2));
yMinIris5 = ceil(n2/2 - n2 * sqrt(0.05 / 2));
yMaxIris5 = ceil(n2/2 + n2 * sqrt(0.05 / 2));

% Construct truncatation windows and multiply them by the original DFT in
% order to truncate.

% 30%.
hLena30 = zeros(m1, n1);
hLena30(xMinLena30:xMaxLena30, yMinLena30:yMaxLena30) = 1;
lenaDFT30 = fftshift(hLena30) .* lenaDFT;

% 15%.
hLena15 = zeros(m1, n1);
hLena15(xMinLena15:xMaxLena15, yMinLena15:yMaxLena15) = 1;
lenaDFT15 = fftshift(hLena15) .* lenaDFT;

% 5%.
hLena5 = zeros(m1, n1);
hLena5(xMinLena5:xMaxLena5, yMinLena5:yMaxLena5) = 1;
lenaDFT5 = fftshift(hLena5) .* lenaDFT;

% 30%.
hIris30 = zeros(m2, n2);
hIris30(xMinIris30:xMaxIris30, yMinIris30:yMaxIris30) = 1;
irisDFT30 = fftshift(hIris30) .* irisDFT;

% 15%.
hIris15 = zeros(m2, n2);
hIris15(xMinIris15:xMaxIris15, yMinIris15:yMaxIris15) = 1;
irisDFT15 = fftshift(hIris15) .* irisDFT;

% 5%.
hIris5 = zeros(m2, n2);
hIris5(xMinIris5:xMaxIris5, yMinIris5:yMaxIris5) = 1;
irisDFT5 = fftshift(hIris5) .* irisDFT;

% Plot truncated DFTs.
fig10 = figure(10);
subplot(3,2,1); imshow(log(abs(fftshift(lenaDFT30))), []); title('Log of 30% Magnitude of Lena Image');
subplot(3,2,3); imshow(log(abs(fftshift(lenaDFT15))), []); title('Log of 15% Magnitude of Lena Image');
subplot(3,2,5); imshow(log(abs(fftshift(lenaDFT5))), []); title('Log of 5% Magnitude of Lena Image');
subplot(3,2,2); imshow(log(abs(fftshift(irisDFT30))), []); title('Log of 30% Magnitude of Iris Image');
subplot(3,2,4); imshow(log(abs(fftshift(irisDFT15))), []); title('Log of 15% Magnitude of Iris Image');
subplot(3,2,6); imshow(log(abs(fftshift(irisDFT5))), []); title('Log of 5% Magnitude of Iris Image');
saveas(fig10, 'fig10_hw6.jpg');

%% d. Apply the 2D Inverse DFT to reconstruct the image for each of the
% truncated spectra. Print out both images reconstructed using each of
% these truncation windows.

% Apply inverse DFT and bring back to original dimensions.
lenaReconstruct = ifft2(lenaDFT);
lenaReconstruct = lenaReconstruct(1:256, 1:256);
lenaReconstruct30 = ifft2(lenaDFT30);
lenaReconstruct30 = lenaReconstruct30(1:256, 1:256);
lenaReconstruct15 = ifft2(lenaDFT15);
lenaReconstruct15 = lenaReconstruct15(1:256, 1:256);
lenaReconstruct5 = ifft2(lenaDFT5);
lenaReconstruct5 = lenaReconstruct5(1:256, 1:256);

irisReconstruct = ifft2(irisDFT);
irisReconstruct = irisReconstruct(1:372, 1:298);
irisReconstruct30 = ifft2(irisDFT30);
irisReconstruct30 = irisReconstruct30(1:372, 1:298);
irisReconstruct15 = ifft2(irisDFT15);
irisReconstruct15 = irisReconstruct15(1:372, 1:298);
irisReconstruct5 = ifft2(irisDFT5);
irisReconstruct5 = irisReconstruct5(1:372, 1:298);

fig11 = figure(11);
subplot(2,2,1); imshow(real(lenaReconstruct), []); title('Reconstructed Image');
subplot(2,2,2); imshow(real(lenaReconstruct30), []); title('Reconstructed 30% Truncated Image');
subplot(2,2,3); imshow(real(lenaReconstruct15), []); title('Reconstructed 15% Truncated Image');
subplot(2,2,4); imshow(real(lenaReconstruct5), []); title('Reconstructed 5% Truncated Image');
saveas(fig11, 'fig11_hw6.jpg');

fig12 = figure(12);
subplot(2,2,1); imshow(real(irisReconstruct), []); title('Reconstructed Image');
subplot(2,2,2); imshow(real(irisReconstruct30), []); title('Reconstructed 30% Truncated Image');
subplot(2,2,3); imshow(real(irisReconstruct15), []); title('Reconstructed 15% Truncated Image');
subplot(2,2,4); imshow(real(irisReconstruct5), []); title('Reconstructed 5% Truncated Image');
saveas(fig12, 'fig12_hw6.jpg');

%% e. Compute SNR for reconstructed images.

lenaSNR30 = real(sum(lenaReconstruct30(:).^2) / sum((lena(:) - lenaReconstruct30(:)).^2));
lenaSNR15 = real(sum(lenaReconstruct15(:).^2) / sum((lena(:) - lenaReconstruct15(:)).^2));
lenaSNR5 = real(sum(lenaReconstruct5(:).^2) / sum((lena(:) - lenaReconstruct5(:)).^2));

irisSNR30 = real(sum(irisReconstruct30(:).^2) / sum((iris(:) - irisReconstruct30(:)).^2));
irisSNR15 = real(sum(irisReconstruct15(:).^2) / sum((iris(:) - irisReconstruct15(:)).^2));
irisSNR5 = real(sum(irisReconstruct5(:).^2) / sum((iris(:) - irisReconstruct5(:)).^2));
