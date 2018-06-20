%% EBME 461 HW07
% Graph Analysis of Pathology Images
% Main Script
% Original Script: @Jacob Antunes

%% Setup
clear all;
close all;
clc;

%% Recurrence
fprintf('\nBiomchemical Recurrence\n');

recurr = imread('0681_A_13_2.tif');

%image
figure(1)
imshow(recurr);
title('biochemical recurrence');

fprintf('\tIdentifying Nuceli...\n');

method = questdlg('How would you like to define the nuclei?','Nuclei Detection','Manual Selection','Automated Segmentation','Automated Segmentation');

if strcmp(method,'Manual Selection')
    [x_recurr, y_recurr] = getpts();
else
    figure;
    recurrbw = rgb2gray(recurr);
    mask = uint8(zeros(size(recurrbw)));
    mask(250:1600,215:1530) = 1;
    mask = mask .* recurrbw;
    mask(find(mask>0 & mask<20)) = 255;
    mask(find(mask<255)) = 0;
    bw = bwmorph(mask,'shrink',inf) + 0;
    imshow(bw == 0);
    title('biochemical recurrence nuclei');
    [y_recurr, x_recurr] = find(bw);
end

%% Nonrecurrence
fprintf('\nNon-recurrence\n');

nonrecurr = imread('0681_A_12_14.tif');

%image
figure;
imshow(nonrecurr);
title('nonrecurrence');

fprintf('\tIdentifying Nuceli...\n');

method = questdlg('How would you like to define the nuclei?','Nuclei Detection','Manual Selection','Automated Segmentation','Automated Segmentation');

if strcmp(method,'Manual Selection')
    [x_nonrecurr, y_nonrecurr] = getpts(figure(5));
else
    figure;
    nonrecurrbw = rgb2gray(nonrecurr);
    mask = uint8(zeros(size(nonrecurrbw)));
    mask(189:1407,267:1530) = 1;
    mask(200:300,400:500) = 0;
    mask = mask .* nonrecurrbw;
    mask(find(mask>15 & mask<52)) = 255;
    mask(find(mask<255)) = 0;
    bw = bwmorph(mask,'shrink',inf) + 0;
    imshow(bw == 0);
    title('biochemical recurrence nuclei');
    [y_nonrecurr, x_nonrecurr] = find(bw);
end

%% Compute Graphs
%Voronoi
fprintf('\tComputing Voronoi for recurrence case...\n');
[vx_recurr, vy_recurr] = voronoi(x_recurr,y_recurr);
figure(2);
imshow(recurr); alpha(0.5);
hold on;
plot(x_recurr,y_recurr,'r+',vx_recurr,vy_recurr,'b-'); hold off;
title('biochemical recurrence: Voronoi');

%Delaunay
fprintf('\tComputing Delaunay for recurrence case...\n');
tri_recurr = delaunay(x_recurr,y_recurr);
figure(3);
imshow(recurr); alpha(0.5);
hold on;
triplot(tri_recurr, x_recurr, y_recurr);
title('biochemical recurrence: Delaunay');


%Voronoi
fprintf('\tComputing Voronoi for nonrecurrence case...\n');
[vx_nonrecurr, vy_nonrecurr] = voronoi(x_nonrecurr, y_nonrecurr);
figure(4);
imshow(nonrecurr); alpha(0.5);
hold on;
plot(x_nonrecurr, y_nonrecurr, 'r+', vx_nonrecurr, vy_nonrecurr, 'b-'); 
hold off;
title('biochemical nonrecurrence: Voronoi');

%Delaunay
fprintf('\tComputing Delaunay for nonrecurrence case...\n');
tri_nonrecurr = delaunay(x_nonrecurr,y_nonrecurr);
figure(5);
imshow(nonrecurr); alpha(0.5);
hold on;
triplot(tri_nonrecurr, x_nonrecurr, y_nonrecurr);
title('biochemical nonrecurrence: Delaunay');
