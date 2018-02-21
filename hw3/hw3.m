%% Data Loading.
close all; clear; clc

load('Atlas.mat');
load('Brain.mat');

%% Control Point Selection.
% Use the cpselect() function to select 8 control points. We will use these
% points to solve for parameters A and B, column vectors containing the
% transformation matrices.

% Use Atlas image as the MOVING and Brain image as the FIXED params.
% [movingPoints, fixedPoints] = cpselect(Atlas, Brain, 'Wait', true);
% save('movingPoints');
% save('fixedPoints');
load('movingPoints.mat');
load('fixedPoints.mat');

%% Matrix Construction.
% Separate (x,y) pairs for the Brain (non-warped, fixed) control points.
X_u = fixedPoints(:,1);
Y_u = fixedPoints(:,2);

% Separate (x,y) pairs for the Atlas (warped, moving) control points.
X_w = movingPoints(:,1);
Y_w = movingPoints(:,2);

% Construct matrix D using (x,y) pairs above.
D = [ones(8,1), X_u, Y_u, X_u .* Y_u];

%% Parameter Estimation.
% Determine transformation parameters using least-square approach.
inverse = pinv(D);
A = inverse * X_w;
B = inverse * Y_w;

fprintf('The parameters of matrix A are:\n%.4f\n%.4f\n%.4f\n%.4f\n', A(1), A(2), A(3), A(4));
fprintf('The parameters of matrix B are:\n%.4f\n%.4f\n%.4f\n%.4f\n', B(1), B(2), B(3), B(4));
%% Affine Transformation.

% Simultaneously produce transformed image and edge image.
AtlasEdges = edge(Atlas, 'canny');

% Retrieve row and column size of the image for looping purposes.
[row, col] = size(Atlas);
newAtlas = zeros(row, col);
newAtlasEdges = zeros(row, col);
for i = 1:row
    for j = 1:col
        % Calculate (x,y) location in the atlas image (xOut, yOut)
        % for the given pixel in the MRI image (i,j).
        % (i, j) is the desired position in the unwarped image.
        x = A(1) + A(2) * i + A(3) * j + A(4) * i * j;
        y = B(1) + B(2) * i + B(3) * j + B(4) * i * j;
        
        % Use the nearest neighbor approach to interpolate grey values.
        x = round(x);
        y = round(y);
        if (x > 0 && y > 0 && x <= row && y <= col)
            newAtlas(j, i) = Atlas(y, x);
            newAtlasEdges(j, i) = AtlasEdges(y, x);
        end
    end
end
newAtlas = uint8(newAtlas);

fig1 = figure(1);
subplot(2,2,1); imshow(Atlas); title('Original Atlas Image');
subplot(2,2,2); imshow(AtlasEdges); title('Original Atlas Edge Image');
subplot(2,2,3); imshow(newAtlas); title('Warped Atlas Image');
subplot(2,2,4); imshow(newAtlasEdges); title('Warped Atlas Edge Image');
saveas(fig1, 'hw3_fig1.jpg')

%% Calculate Euclidean Distance between Control Points
% Calculate the estimated point pair.
X_w_calc = D * A;
Y_w_calc = D * B;

sse_x = sum((X_w_calc - X_u).^2);
sse_y = sum((Y_w_calc - Y_u).^2);
avgDist = sqrt(sse_x + sse_y) / 8;
fprintf('The average Euclidean distance is %.4f\n', avgDist);

%% Create Edge Image of Atlas Image Superimposed on transformed MRI image.
overlayImage = imoverlay(Brain, newAtlasEdges);
fig2 = figure(2);
imshow(overlayImage); title('Contoured Brain Image')
saveas(fig2, 'hw3_fig2.jpg');

%% Extra Credit: Create a Superfused image.

% Set the factor to 0.5 to mix the two images equally.
alphaFactor = 0.6;

bgImg = double(Brain);
fgImg = double(newAtlas);

bgImageAlpha = (1 - alphaFactor) .* bgImg;
fgImageAlpha = alphaFactor .* fgImg;

fusedImg = bgImageAlpha + fgImageAlpha;

fig3 = figure(3);
imshow(uint8(fusedImg)); title('Superfused Image');
saveas(fig3, 'hw3_fig3.jpg');

%% Extra Credit: Create an Edge Image to Superimpose (green/blue).