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
%% 4.2.1, 4.2.2
img2 = imread("./assets/Image02.jpg");
noisy = imnoise(img2, 'Gaussian', 0, .04);

figure('name', "Raw Vs. Noisy");
subplot(1,2,1)
imshow(img2);
title('Original Image');
subplot(1,2,2)
imshow(noisy);
title('Noisy Image');

%% 4.2.3
kernel = ones(3,3) / 9;
filtered = imfilter(noisy,kernel);

figure('name', "Filtered Image");
subplot(1,3,1)
imshow(img2);
title('Original Image');
subplot(1,3,2)
imshow(noisy);
title('Noisy Image');
subplot(1,3,3)
imshow(filtered);
title('Filtered Image');

%% 4.2.4
kernel = ones(5,5) / 25;
filtered = imfilter(noisy,kernel);

figure('name', "Filtered Image");
subplot(1,3,1)
imshow(img2);
title('Original Image');
subplot(1,3,2)
imshow(noisy);
title('Noisy Image');
subplot(1,3,3)
imshow(filtered);
title('Filtered Image');
%% 4.2.5
salt_papper_noisy = imnoise(img2, 'salt & pepper', 0.1);

figure('name', "Raw Vs. Noisy");
subplot(1,2,1)
imshow(img2);
title('Original Image');
subplot(1,2,2)
imshow(salt_papper_noisy);
title('Noisy Image');
