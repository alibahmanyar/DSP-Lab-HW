clc;
clear all;
close all;
%% 4.1.1
img = imread("./assets/lena.png");

imshow(img)

%% 4.1.2
img = im2double(img);

%% 4.1.3
figure('name', "Lenas Histogram")
imhist(img)
%% 4.1.4
enhanced = histeq(img);
figure('name', "Enhanced Contrast");
subplot(1,2,1)
imshow(img(:,:,:));
title('Slice of Original Image');
subplot(1,2,2)
imshow(enhanced(:,:,:));
title('Slice of Enhanced Image');

figure('name', "montage pair");
imshowpair(img,enhanced,'montage');
%% 4.1.5
figure('name', "Enhanced Lena Histogram")
imhist(enhanced)
%%
