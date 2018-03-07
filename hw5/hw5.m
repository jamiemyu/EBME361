%% Problem 1

A = [0 1 0; 0 1 0; 0 1 0];
B = [0 0 0; 1 1 1; 0 0 0];
C = [0 1 0; 1 1 1; 0 1 0];

output1 = imdilate(A, B, 'full');
disp(output1);

output2 = imdilate(output1, C, 'full');
disp(output2);

%% Problem 2
close all; clear; clc;

% Load the image.
cameraMan = imread('cameraman.png');

fig1 = figure(1);
imshow(cameraMan); title('Original Image');

% a. Create a structuring element to dilate image to remove pepper noise.
SE_dilate = strel('rectangle', [3, 2]);
dilatedCameraMan = imdilate(cameraMan, SE_dilate);
fig2 = figure(2);
imshow(dilatedCameraMan); title('Dilated Image');

% b. Design an open-close gray scale morphology operation to reduce both
% salt and pepper noise. We perform consecutive closing and opening
% operations on the image to first remove the pepper (close the "holes" in
% the image and then the salt (open the "noise").
SE = strel('rectangle', [4, 4]);
closedCameraMan = imclose(cameraMan, SE);
finalCameraMan = imopen(closedCameraMan, SE);
fig3 = figure(3);
imshow(finalCameraMan); title('Salt/Pepper Removed Image');

saveas(fig1, 'hw5_fig1.jpg');
saveas(fig2, 'hw5_fig2.jpg');
saveas(fig3, 'hw5_fig3.jpg');

%% Problem 3
close all; clear; clc;

% Load the image.
circles = imread('Circles.png');

fig4 = figure(4);
imshow(circles); title('Original Image');
saveas(fig4, 'hw5_fig4.jpg');

% Perform thresholding of the image. First convert rgb image to gray, then
% perform gray thresholding.
grayCircles = rgb2gray(circles);

% a. Define threshold at pixel value 90 to create a binary image.
threshold = 90;
thresholdCircles = grayCircles < threshold;

fig5 = figure(5);
imshow(thresholdCircles); title('Threshold Image');
saveas(fig5, 'hw5_fig5.jpg');

% b. Use erosion operation to separate the circles distinctly.
SE_erode = strel('disk', 11);
erodedCircles = imerode(thresholdCircles, SE_erode);

fig6 = figure(6);
imshow(erodedCircles); title('Separated Circles');
saveas(fig6, 'hw5_fig6.jpg');

%% Problem 4
 clear; clc;

% Load the image.
circlesLines = imread('Circle_and_Lines.png');

% Convert to 2D image (gray-scale);
circlesLines = rgb2gray(circlesLines);

% Find the threshold.
circleLevel = graythresh(circlesLines);

% Binarize the image.
circlesLines = imbinarize(circlesLines, circleLevel);

fig7 = figure(7);
imshow(circlesLines); title('Original Image');
saveas(fig7, 'hw5_fig7.jpg');

circlesImage = imopen(circlesLines, strel('disk', 19));
lines1 = imopen(circlesLines, strel('line', 150, 0));
lines2 = imopen(circlesLines, strel('line', 145, 275));
lines3 = imopen(circlesLines, strel('line', 200, 20));
lines4 = imopen(circlesLines, strel('line', 200, 77));
lines5 = imopen(circlesLines, strel('line', 110, 135));
lines6 = imopen(circlesLines, strel('line', 160, -15));
lines7 = imopen(circlesLines, strel('line', 200, 5));

horizontalLines = lines1 | lines5 | lines6 | lines7;
verticalLines = lines2 | lines3 | lines4;

linesImage = horizontalLines | verticalLines;
linesImage = imopen(linesImage, strel('disk', 7));

fig8 = figure(8);
subplot(1,2,1); imshow(circlesImage); 
title('Circles Image');
subplot(1,2,2); imshow(linesImage); 
title('Lines Image');
saveas(fig8, 'hw5_fig8.jpg');

%% c. Compare different SE sizes.
SE_open = strel('disk', 13);
circlesImage = imopen(circlesLines, SE_open);

fig9 = figure(9);
subplot(2,2,1); imshow(circlesImage);
title('Circles Image opened with size 13');

SE_open = strel('disk', 16);
circlesImage = imopen(circlesLines, SE_open);

fig9;
subplot(2,2,2); imshow(circlesImage);
title('Circles Image opened with size 16');

SE_open = strel('disk', 27);
circlesImage = imopen(circlesLines, SE_open);

fig9;
subplot(2,2,3); imshow(circlesImage);
title('Circles Image opened with size 27');

SE_open = strel('disk', 32);
circlesImage = imopen(circlesLines, SE_open);

fig9;
subplot(2,2,4); imshow(circlesImage);
title('Circles Image opened with size 32');

saveas(fig9, 'hw5_fig9.jpg');

%% d. Write an algorithm to count the number of items in an image.
close all; clear; clc; 

circlesLines = imread('Circle_and_Lines.png');
circlesLines = rgb2gray(circlesLines);
falconLevel = graythresh(circlesLines);
circlesLines = imbinarize(circlesLines, falconLevel);

SE_open = strel('disk', 19);
circlesImage = imopen(circlesLines, SE_open);
lines1 = imopen(circlesLines, strel('line', 150, 0));
lines2 = imopen(circlesLines, strel('line', 145, 275));
lines3 = imopen(circlesLines, strel('line', 200, 20));
lines4 = imopen(circlesLines, strel('line', 200, 77));
lines5 = imopen(circlesLines, strel('line', 110, 135));
lines6 = imopen(circlesLines, strel('line', 160, -15));
lines7 = imopen(circlesLines, strel('line', 200, 5));

horizontalLines = lines1 | lines5 | lines6 | lines7;
verticalLines = lines2 | lines3 | lines4;

linesImage = horizontalLines | verticalLines;
linesImage = imopen(linesImage, strel('disk', 7));

% Label the elements in each image.
circleLabels = logical(circlesImage);

s1 = regionprops(circleLabels, 'Centroid');
fig10 = figure(10);
subplot(1,2,1); imshow(circlesImage); title('Labeled Circles Image');
hold on
for j = 1:numel(s1)
    c = s1(j).Centroid;
    text(c(1), c(2), sprintf('%d', j), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end
hold off

fprintf('There are %d circles in the image.\n', numel(s1));
saveas(fig10, 'hw5_fig10.jpg');

% Label the elements in each image.
verticalLinesLabels = logical(verticalLines);
horizontalLinesLabels = logical(horizontalLines);

s2 = regionprops(horizontalLinesLabels, 'Centroid');
fig10;
subplot(1,2,2); imshow(linesImage); title('Labeled Lines Image');
hold on
for j = 1:numel(s2)
    c = s2(j).Centroid;
    text(c(1), c(2), sprintf('%d', j), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 5);
end

s3 = regionprops(verticalLinesLabels, 'Centroid');
for k = 1:numel(s3)
    i = k + j;
    c = s3(k).Centroid;
    text(c(1), c(2), sprintf('%d', i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'FontSize', 5);
end
hold off

fprintf('There are %d lines in the image.\n', numel(s2) + numel(s3));
saveas(fig10, 'hw5_fig10.jpg');
%% Problem 5
close all; clear; clc;

% Load images and perform thresholding.
falcon = imread('Falcon.png');
falcon = rgb2gray(falcon);
falconLevel = graythresh(falcon);

girl = imread('Girl.png');
girl = rgb2gray(girl);
girlLevel = graythresh(girl);

dog = imread('dog.png');
dog = rgb2gray(dog);
dogLevel = graythresh(dog);

balloon = imread('Balloon.png');
balloon = rgb2gray(balloon);
balloonLevel = graythresh(balloon);

% Binarize the images.
falcon = imcomplement(imbinarize(falcon, falconLevel));
girl = imcomplement(imbinarize(girl, girlLevel));
dog = imcomplement(imbinarize(dog, dogLevel));
balloon = imcomplement(imbinarize(balloon, balloonLevel));

% Zero-pad the left side so the balloon is not touching the edges of the
% image.
[r, ~] = size(balloon);
balloon = [false(r,20), balloon];
[~, c] = size(balloon);
balloon = [false(20,c); balloon];

fig14 = figure(14);
subplot(2,2,1); imshow(falcon); title('Original Falcon Image');
subplot(2,2,2); imshow(girl); title('Original Girl Image');
subplot(2,2,3); imshow(dog); title('Original Dog Image');
subplot(2,2,4); imshow(balloon); title('Original Balloon Image');
saveas(fig14, 'hw5_fig14.jpg');

% Perform closing of images with disk shaped structuring element.
SE_close_med = strel('disk', 10);
SE_close_large = strel('disk', 20);

falconClosed = imclose(falcon, SE_close_med);
girlClosed = imclose(girl, SE_close_large);
dogClosed = imclose(dog, SE_close_med);
balloonClosed = imclose(balloon, SE_close_med);

fig15 = figure(15);
subplot(2,2,1); imshow(falconClosed); title('Closed Falcon Image with SE size 10');
subplot(2,2,2); imshow(girlClosed); title('Closed Girl Image with SE size 20');
subplot(2,2,3); imshow(dogClosed); title('Closed Dog Image with SE size 10');
subplot(2,2,4); imshow(balloonClosed); title('Closed Balloon Image with SE size 10');
saveas(fig15, 'hw5_fig15.jpg');

%% a. Create skeleton images from closed images.
close all; 

skeletonizedFalcon = bwmorph(falconClosed, 'skel', Inf);
skeletonizedGirl = bwmorph(girlClosed, 'skel', Inf);
skeletonizedDog = bwmorph(dogClosed, 'skel', Inf);
skeletonizedBalloon = bwmorph(balloonClosed, 'skel', Inf);

fig16 = figure(16); 
imshow(skeletonizedFalcon); title('Skeleton Falcon Image');
fig17 = figure(17); 
imshow(skeletonizedGirl); title('Skeleton Girl Image');
fig18 = figure(18); 
imshow(skeletonizedDog); title('Skeleton Dog Image');
fig19 = figure(19); 
imshow(skeletonizedBalloon); title('Skeleton Balloon Image');

saveas(fig16, 'hw5_fig16.jpg');
saveas(fig17, 'hw5_fig17.jpg');
saveas(fig18, 'hw5_fig18.jpg');
saveas(fig19, 'hw5_fig19.jpg');

%% b. Compute the center of mass of the shapes of the four images.
close all;

% Use shrinking operation with infinite iterations.
comFalcon = bwmorph(skeletonizedFalcon, 'shrink', Inf);
comGirl = bwmorph(skeletonizedGirl, 'shrink', Inf);
comDog = bwmorph(skeletonizedDog, 'shrink', Inf);
comBalloon = bwmorph(skeletonizedBalloon, 'shrink', Inf);

% Find the indices of center of mass locations.
[mFalcon, nFalcon] = find(comFalcon == 1);
[mGirl, nGirl] = find(comGirl == 1);
[mDog, nDog] = find(comDog == 1);
[mBalloon, nBalloon] = find(comBalloon == 1);

% Plot the center of mass on each of the images.
fig20 = figure(20); subplot(2,2,1);
imshow(falcon); hold on; plot(nFalcon, mFalcon, 'r.', 'MarkerSize', 15);
fig20; subplot(2,2,2);
imshow(girl); hold on; plot(nGirl, mGirl, 'r.', 'MarkerSize', 15);
fig20; subplot(2,2,3);
imshow(dog); hold on; plot(nDog, mDog, 'r.', 'MarkerSize', 15);
fig20; subplot(2,2,4);
imshow(balloon); hold on; plot(nBalloon, mBalloon, 'r.', 'MarkerSize', 15);
saveas(fig20, 'hw5_fig20.jpg');

%% c. Compute the area of each of the shape of images.
close all;

areaFalcon = bwarea(falconClosed);
disp(areaFalcon);
areaGirl = bwarea(girlClosed);
disp(areaGirl);
areaDog = bwarea(dogClosed);
disp(areaDog);
areaBalloon = bwarea(balloonClosed);
disp(areaBalloon);
