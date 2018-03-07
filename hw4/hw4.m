%% Load matrices into the workspace. 

load('Contrast1.mat');
load('Contrast2.mat');

%% Problem 1
% Write function myNCC(img1, img2) which computes normalized cross
% correlation between two input images.

%% Problem 2
close all; clc;

% Define translation parameters.
tx = -1.8;
ty = -2.1;

X = 1:length(Contrast2);
Y = 1:length(Contrast2);
[X, Y] = meshgrid(X, Y);

Xi = X + tx;
Yi = Y + ty;
translatedContrast2 = interp2(X, Y, Contrast2, Xi, Yi, 'linear', 0);

fig1 = figure(1); subplot(1,2,1)
imagesc(translatedContrast2); 
colormap gray; axis image; axis off; title('Translated Contrast2 Image');

% Create a subtracted image to confirm registration.
% This image should be black if registration was correct.
subtractedImage = Contrast1 - translatedContrast2;
fig1; subplot(1,2,2)
imagesc(subtractedImage);
colormap gray; axis image; axis off; title('Subtracted Image');
saveas(fig1, 'hw4_fig1.jpg');

% For the purpose of computing NCC correctly, set the edge values that had 
% no points mapped to them equal to the nearest valued row and column. 
% Otherwise, the NCC calculation will be skewed since the correlation 
% between the rows andcolumns of zeros will be distinctly different than 
% the same rows and columns of the fixed image.
translatedContrast2(1:3,:) = [translatedContrast2(4,:); ...
                              translatedContrast2(4,:);
                              translatedContrast2(4,:)];
translatedContrast2(:,1:2) = [translatedContrast2(:,3), ...
                              translatedContrast2(:,3)];
                          
% Compute NCC.
R = myNCC(Contrast1, translatedContrast2);
fprintf('The normalized cross correlation between Contrast1 and registered Contrast2 is %.4f\n', R);

%% Problem 3
% Register the same image pair using fminsearch.
close all; clc;

% Start initial transformation parameters at tx = 0, ty = 0 if only 
% translation is used.
p0 = [0 0];

% Define options e.g. (doc optimset for details).
opt = optimset('Display', 'iter');

% Actual call to optimization.
scalingParam = 1;
[phat, fval] = fminsearch(@funRegister, p0, opt, Contrast1, Contrast2, scalingParam);

% The cross correlation sign is switched again since we calculated it for a
% minimum in the fminsearch() function
nccFound = -fval;
fprintf('The normalized correlation coefficient was found to be %.4f\n', nccFound);

txFound = phat(1);
tyFound = phat(2);

% Compute xi and yi from x, y, and phat.
n = length(Contrast2);
[X, Y] = meshgrid(1:n, 1:n);
xi = X + txFound;
yi = Y + tyFound;

newContrast2 = interp2(X, Y, Contrast2, xi, yi, 'linear', 0);
subtractedImage = Contrast1 - newContrast2;
fig2 = figure(2);
imagesc(subtractedImage);
colormap gray; axis image; axis off; title('Subtracted Image');
saveas(fig2, 'hw4_fig2.jpg');

fig3 = figure(3);
subplot(1,2,1);
imagesc(Contrast1);
colormap gray; axis image; axis off; title('Contrast1 Image');
subplot(1,2,2);
imagesc(newContrast2); 
colormap gray; axis image; axis off; title('Translated Contrast2 Image');
saveas(fig3, 'hw4_fig3.jpg');

%% Problem 4
% Register the same image pair using fminsearch and modified parameters.
%close all; clc;

% Start initial transformation parameters at tx = 0, ty = 0 if only 
% translation is used.
p0 = [0 0];

% Define options e.g. (doc optimset for details)
opt = optimset('Display', 'iter', 'PlotFcns',@optimplotfval);

% Actual call to optimization. Pass the scaling parameter of choice.
scalingParam = 1;
[phat, fval] = fminsearch(@funRegister, p0, opt, Contrast1, Contrast2, scalingParam);

% The cross correlation sign is switched again since we calculated it for a
% minimum in the fminsearch() function
nccFound = -fval;
fprintf('The normalized correlation coefficient was found to be %.4f\n', nccFound);

% Multiply by the scaling parameter.
txFound = scalingParam * phat(1);
tyFound = scalingParam * phat(2);

% Compute xi and yi from x, y, and phat.
n = length(Contrast2);
[X, Y] = meshgrid(1:n, 1:n);
xi = X + txFound;
yi = Y + tyFound;

newContrast2 = interp2(X, Y, Contrast2, xi, yi, 'linear', 0);

fig4 = figure(4); 
imagesc(Contrast1 - newContrast2);
colormap gray; axis image; axis off; title('Subtract Image');
saveas(fig4, 'hw4_fig4.jpg');

%% Problem 5
% Register a new image pair using fminsearch and modified parameters.
close all; clear; clc;

load('mri1.mat');
load('mri2.mat');

fig5 = figure(5);
subplot(1,2,1); imagesc(mri1 - mri2);
colormap gray; axis image; axis off; title('Subtract Image of MRI');

% Start initial transformation parameters at tx = 0, ty = 0 if only 
% translation is used.
p0 = [0 0];

% Define options e.g. (doc optimset for details)
opt = optimset('Display', 'iter', 'TolX', 1, 'TolFun', 1e-3, 'PlotFcns',@optimplotfval);

% Actual call to optimization. Pass the scaling parameter of choice.
scalingParam = 1000;
[phat, fval] = fminsearch(@funRegister, p0, opt, mri1, mri2, scalingParam);

% The cross correlation sign is switched again since we calculated it for a
% minimum in the fminsearch() function
nccFound = -fval;
fprintf('The normalized correlation coefficient was found to be %.4f\n', nccFound);

% Multiply by the scaling parameter.
txFound = scalingParam * phat(1);
tyFound = scalingParam * phat(2);

% Compute xi and yi from x, y, and phat.
n = length(mri2);
[X, Y] = meshgrid(1:n, 1:n);
xi = X + txFound;
yi = Y + tyFound;

newMri2 = interp2(X, Y, mri2, xi, yi, 'linear', 0);

figure(5); 
subplot(1,2,2); imagesc(mri1 - newMri2);
colormap gray; axis image; axis off; title('Subtract Image of MRI');
saveas(fig5, 'hw4_fig5.jpg');

fig6 = figure(6);
subplot(1,2,1); imagesc(mri1);
colormap gray; axis image; axis off; title('mri1');
subplot(1,2,2); imagesc(newMri2);
colormap gray; axis image; axis off; title('Translated mri2');
saveas(fig6, 'hw4_fig6.jpg');

%% Problem 6 (Extra Credit)
% Register a new image pair using fminsearch and modified parameters.
close all; clear; clc;

load('origmri.mat');
load('taskmri.mat');

fig7 = figure(7);
subplot(1,2,1); imagesc(origmri - taskmri);
colormap gray; axis image; axis off; title('Subtract Image of MRI before registration');

% Start initial transformation parameters at tx = 0, ty = 0 if only 
% translation is used.
p0 = [0 0 10];

% Define options e.g. (doc optimset for details)
opt = optimset('Display', 'iter', 'TolX', 1e-2, 'TolFun', 1e-3, 'PlotFcns', @optimplotfval);

% Actual call to optimization. Pass the scaling parameter of choice.
scalingParam = 100;
[phat, fval] = fminsearch(@funRegisterRot, p0, opt, origmri, taskmri, scalingParam);

% The cross correlation sign is switched again since we calculated it for a
% minimum in the fminsearch() function
nccFound = -fval;
fprintf('The normalized correlation coefficient was found to be %.4f\n', nccFound);

% Multiply by the scaling parameter.
txFound = scalingParam * phat(1);
tyFound = scalingParam * phat(2);
thetaFound = scalingParam * phat(3);

% Compute xi and yi from x, y, and phat.
n = length(taskmri);
[X, Y] = meshgrid(1:n, 1:n);
xi = X + txFound;
yi = Y + tyFound;

newTaskMri = interp2(X, Y, taskmri, xi, yi, 'linear', 0);

% Rotate image.
newTaskMri = imrotate(newTaskMri, deg2rad(thetaFound), 'bilinear', 'crop');

figure(7); 
subplot(1,2,2); imagesc(origmri - newTaskMri);
colormap gray; axis image; axis off; title('Subtract Image of MRI after registration');
saveas(fig7, 'hw4_fig7.jpg');

fig8 = figure(8);
subplot(1,2,1); imagesc(origmri);
colormap gray; axis image; axis off; title('Original MRI');
subplot(1,2,2); imagesc(newTaskMri);
colormap gray; axis image; axis off; title('Registered Task MRI');
saveas(fig8, 'hw4_fig8.jpg');