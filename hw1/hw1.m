%% Problem 2: Unsharp Mask Edge Enhancement
close all; clear; clc;

% Load the data.
f = imread('lena.tif');
fOrg = imread('lena_org.tif');

% Show the image.
fig1 = figure(1); subplot(1,3,1); imshow(f); title('lena degraded image');
saveas(fig1, 'problem2_fig1.jpg');
figure(1); subplot(1,3,2); imshow(fOrg); title('original lena image');

% Construct a 5x5 filter kernel.
w = fspecial('average', 5);

% Convolve the input with the smoothing window.
fConv = imfilter(f, w, 'conv');

% Select a beta value.
beta = 0.5;
% Construct the new image using the difference image.
g = f + beta * (f - fConv);
figure(1); subplot(1,3,3); imshow(g); title('sharp output image');
saveas(fig1, 'problem2_fig1.jpg');

%% Problem 3: Sobel edge filtering using imfilter()
clear; close all; clc;
% Load the data.
f = imread('TUMOR.tif');
fig3 = figure(3); imshow(f); title('original tumor image');

% Apply the Sobel filter to the image.
outputImg = edge(f(:,:,1), 'sobel', 0.13);
fig4 = figure(4); imshow(outputImg); title('Sobel filtered image');
saveas(fig4, 'problem3_fig1.jpg');

%% Problem 4: Investigate edge detection on the image, EDGES.mat
close all; clear; clc;
f = imread('EDGES.tif');
figure(5); imshow(f);

% Sobel
[sobelDisk, sobelThresh] = edge(f(:,:,1), 'Sobel');
fig6 = figure(6); imshow(sobelDisk); title('Sobel');
saveas(fig6,'problem4_fig1.jpg');

% Prewitt
[prewittDetect, prewittThresh] = edge(f(:,:,1), 'Prewitt');
fig7 = figure(7); imshow(prewittDetect); title('Prewitt');
saveas(fig7,'problem4_fig2.jpg');

% Roberts
[robertsDetect, robertsThresh] = edge(f(:,:,1), 'Roberts');
fig8 = figure(8); imshow(robertsDetect); title('Roberts');
saveas(fig8,'problem4_fig3.jpg');

% Laplacian of Gaussian (log)
[logDetect, logThresh] = edge(f(:,:,1), 'log');
fig9 = figure(9); imshow(logDetect); title('log');
saveas(fig9,'problem4_fig4.jpg');

% Zero crossing
[zeroCrossDetect, zeroCrossThresh] = edge(f(:,:,1), 'zerocross');
fig10 = figure(10); imshow(zeroCrossDetect); title('zero-cross');
saveas(fig10,'problem4_fig10.jpg');

% Canny
[cannyDetect, cannyThresh] = edge(f(:,:,1), 'Canny');
fig11 = figure(11); imshow(cannyDetect); title('Canny');
saveas(fig11,'problem4_fig11.jpg');

%% Problem 4, Part c
close all; clear; clc;
disk_img = imread('disk_img.tif');
fig15 = figure(15); imshow(disk_img);
saveas(fig15, 'problem4_fig7.jpg');
fig17 = figure(17); imshow(permute(disk_img);

% Sobel
[sobelDisk, sobelThresh] = edge(disk_img(:,:,1), 'Sobel');
fig16 = figure(16); imshow(sobelDisk); title('Sobel Disk');
saveas(fig16, 'problem4_fig8.jpg');
%% Problem 7: Averaging Filter Obeying Linearity
close all; clear; clc;

% Construct a 7x7 averaging filter.
h = fspecial('average', 7);

% Construct a sinusoid test image.
Fs = 0.005; 
dim = 1000; 

I1 = zeros(dim);
for m = 0 : dim - 1 
    for n = 0 : dim - 1
        I1(m + 1, n + 1) = sin(2 * pi * Fs * m);
    end
end
I1 = I1';

% Show the sinusoid.
colormap(gray(256));
fig12 = figure(12); imshow(I1); title('I1');
saveas(fig12, 'problem7_fig1.jpg');

% Construct an image with rectangles and disks.
I2 = zeros(dim);

for i = 200:400
    for j = 200:450
        I2(i,j) = 1;
    end
end

for i = 600:800
    for j = 500:650
        I2(i,j) = 1;
    end
end

imageSizeX = 320;
imageSizeY = 240;
[columnsInImage rowsInImage] = meshgrid(1:324, 1:240);
centerX = 160;
centerY = 120;
radius = 50;
circlePixels = (rowsInImage - centerY).^2 ...
    + (columnsInImage - centerX).^2 <= radius.^2;
circlePixels = double(circlePixels);

I2(580:819, 70:393) = circlePixels;
I2(70:309, 580:903) = circlePixels;

colormap(gray(256));
fig13 = figure(13); imshow(I2); title('I2');
saveas(fig13, 'problem7_fig2.jpg');

% Define scaling parameter (arbitrary).
a = 10;
b = 5;

% Perform f(a*I1 + b*I2) (left-side of the equation).
innerSum = a * I1 + b * I2;
IConv1 = imfilter(innerSum, h, 'conv');

% Perform a*f(I1) + b*f(I2) (right-side of the equation).
IConv2 = a * imfilter(I1, h, 'conv') + b * imfilter(I2, h, 'conv');

% Compare the two filtered images.
fig14 = figure(14); 
subplot(1,2,1); imshow(IConv1); title('left-side equation');
subplot(1,2,2); imshow(IConv2); title('right-side equation');
saveas(fig14, 'problem7_fig3.jpg');