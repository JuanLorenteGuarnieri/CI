%%%%%%%%%%%%%%%%%%%     2.1 Reading the image into Matlab     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the image in TIFF format
workspacePath = fileparts(mfilename('fullpath'));
imgName = 'IMG_1026';
imgPath = fullfile(workspacePath, '/images_tiff/', [imgName, '.tiff']);
cpu_img = imread(imgPath);
img = gpuArray(cpu_img);

% Check and report how many bits per pixel the image has
info = imfinfo(imgPath);
bitsPerPixel = info.BitDepth;
fprintf('Bits per pixel: %d\n', bitsPerPixel);

% Check its width and its height
[height, width, ~] = size(img);
fprintf('Width: %d, Height: %d\n', width, height);

% Convert the image into a double-precision array
img_double = double(img);

% Display the original image
% figure;
% imshow(img);
% title('Original Image');

% Display the double-precision image
% figure;
% imshow(img_double, []);
% title('Double-Precision Image');

%%%%%%%%%%%%%%%%%%%             2.2 Linearization             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linearize the image
black_level = 1023;
saturation_level = 15600;

% Apply linear transformation
img_linear = (img_double - black_level) / (saturation_level - black_level);

% Clip values to the range [0, 1]
img_linear = max(0, min(1, img_linear));

% Display the linearized image
% figure;
% imshow(img_linear, []);
% title('Linearized Image');

%%%%%%%%%%%%%%%%%%%               3 Demosaicing               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify the Bayer pattern
bayer_pattern = 'grbg'; % Replace with the identified pattern

% Separate the color channels based on the Bayer pattern
red_channel = img_linear(1:2:end, 2:2:end);
green_channel1 = img_linear(1:2:end, 1:2:end);
green_channel2 = img_linear(2:2:end, 2:2:end);
blue_channel = img_linear(2:2:end, 1:2:end);

% Interpolate missing pixels using nearest neighbor
red_channel_full = imresize(red_channel, 2, 'nearest');
green_channel_full = imresize((green_channel1 + green_channel2) / 2, 2, 'nearest');
blue_channel_full = imresize(blue_channel, 2, 'nearest');

% Combine the channels to form the demosaiced image
img_demosaiced_nn = cat(3, red_channel_full, green_channel_full, blue_channel_full);

% Display the demosaiced image using nearest neighbor interpolation
% figure;
% imshow(img_demosaiced_nn, []);
% title('Demosaiced Image (Nearest Neighbor)');

% Interpolate missing pixels using bilinear interpolation
red_channel_full = imresize(red_channel, 2, 'bilinear');
green_channel_full = imresize((green_channel1 + green_channel2) / 2, 2, 'bilinear');
blue_channel_full = imresize(blue_channel, 2, 'bilinear');

% Combine the channels to form the demosaiced image
img_demosaiced_bilinear = cat(3, red_channel_full, green_channel_full, blue_channel_full);

% Display the demosaiced image using bilinear interpolation
% figure;
% imshow(img_demosaiced_bilinear, []);
% title('Demosaiced Image (Bilinear Interpolation)');

%%%%%%%%%%%%%%%%%%%             4 White balancing             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% White World White Balancing
maxR = max(max(img_demosaiced_bilinear(:,:,1)));
maxG = max(max(img_demosaiced_bilinear(:,:,2)));
maxB = max(max(img_demosaiced_bilinear(:,:,3)));

img_white_world = img_demosaiced_bilinear;
img_white_world(:,:,1) = img_white_world(:,:,1) / maxR;
img_white_world(:,:,2) = img_white_world(:,:,2) / maxG;
img_white_world(:,:,3) = img_white_world(:,:,3) / maxB;

% Display the white world white balanced image
% figure;
% imshow(img_white_world, []);
% title('White World White Balanced Image');

% Gray World White Balancing
meanR = mean(mean(img_demosaiced_bilinear(:,:,1)));
meanG = mean(mean(img_demosaiced_bilinear(:,:,2)));
meanB = mean(mean(img_demosaiced_bilinear(:,:,3)));

img_gray_world = img_demosaiced_bilinear;
img_gray_world(:,:,1) = img_gray_world(:,:,1) / meanR;
img_gray_world(:,:,2) = img_gray_world(:,:,2) / meanG;
img_gray_world(:,:,3) = img_gray_world(:,:,3) / meanB;

% Display the gray world white balanced image
% figure;
% imshow(img_gray_world, []);
% title('Gray World White Balanced Image');

% Manual White Balancing
% Assume the neutral object is at (x, y) in the image
x = 465; % Replace with actual x coordinate
y = 2525; % Replace with actual y coordinate

R = img_demosaiced_bilinear(y, x, 1);
G = img_demosaiced_bilinear(y, x, 2);
B = img_demosaiced_bilinear(y, x, 3);

SR = (R + G + B) / (3 * R);
SG = (R + G + B) / (3 * G);
SB = (R + G + B) / (3 * B);

img_manual = img_demosaiced_bilinear;
img_manual(:,:,1) = img_manual(:,:,1) * SR;
img_manual(:,:,2) = img_manual(:,:,2) * SG;
img_manual(:,:,3) = img_manual(:,:,3) * SB;

% Display the manually white balanced image
figure;
imshow(img_manual, []);
title('Manual White Balanced Image');

%%%%%%%%%%%%%%%%%%%                5 Denoising                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Define Gaussian kernel
sigma = 2;
kernel_size = 6 * sigma + 1;
x = -floor(kernel_size / 2):floor(kernel_size / 2);
gaussian_kernel = exp(-x.^2 / (2 * sigma^2));
gaussian_kernel = gaussian_kernel / sum(gaussian_kernel);

% Apply Gaussian smoothing
img_gaussian = zeros(size(img_manual));
for i = 1:3
  img_gaussian(:,:,i) = conv2(img_manual(:,:,i), gaussian_kernel, 'same');
  img_gaussian(:,:,i) = conv2(img_gaussian(:,:,i), gaussian_kernel', 'same');
end

% % Define median filter kernel size
% median_kernel_size = 3;

% % Apply median filtering using custom function
% cpu_img_manual = gather(img_manual);
% img_median = zeros(size(cpu_img_manual));
% for i = 1:3
%   img_median(:,:,i) = customMedianFilter(cpu_img_manual(:,:,i), median_kernel_size);
% end
% img_median = gpuArray(img_median);

% Define mean filter kernel
mean_kernel = ones(median_kernel_size) / (median_kernel_size^2);

% Apply mean filtering
img_mean = zeros(size(img_manual));
for i = 1:3
  img_mean(:,:,i) = conv2(img_manual(:,:,i), mean_kernel, 'same');
end


%%%%%%%%%%%%%%%%%%%              6 Color balance              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert the image to HSV color space
img_hsv = rgb2hsv(img_gaussian);

% Boost the saturation channel
saturation_boost = 1.25; % Adjust this value as needed
img_hsv(:,:,2) = img_hsv(:,:,2) * saturation_boost;

% Clip the saturation values to the range [0, 1]
img_hsv(:,:,2) = min(1, img_hsv(:,:,2));

% Convert the image back to RGB color space
img_color_balanced = hsv2rgb(img_hsv);

% Display the color balanced image
figure;
imshow(img_color_balanced, []);
title('Color Balanced Image');

%%%%%%%%%%%%%%%%%%%            7 Tone reproduction            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get gray color
gray = rgb2gray(img_color_balanced);

% Get max grey value
max_gray = max(gray(:));

% Get percentage of the grey
grey_percentage = 0.15;
img_tone_reproduction = img_color_balanced + max_gray*grey_percentage;

% Apply gamma
gamma = 1.75;
img_tone_reproduction(img_tone_reproduction < 0.0031308) = 12.92 * img_tone_reproduction(img_tone_reproduction < 0.0031308);
img_tone_reproduction(img_tone_reproduction >= 0.0031308) = 1.055 * img_tone_reproduction(img_tone_reproduction >= 0.0031308).^(1 / gamma) - 0.055;

% Display the tone reproduced image
figure;
imshow(img_tone_reproduction, []);
title('Tone Reproduced Image');
%%%%%%%%%%%%%%%%%%%               8 Compression               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Export to png
imgPath = fullfile(workspacePath, '/images_output/', [imgName, '.png']);
imwrite(img_tone_reproduction, imgPath);

% Export to jpeg with compression
imgPath = fullfile(workspacePath, '/images_output/', [imgName, '.jpg']);
imwrite(img_tone_reproduction, imgPath, 'Quality', 10);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Custom median filter function
function output = customMedianFilter(input, kernel_size)
pad_size = floor(kernel_size / 2);

% Add padding
padded_input = zeros(size(input) + 2 * pad_size);
padded_input(pad_size+1:end-pad_size, pad_size+1:end-pad_size) = input;

output = zeros(size(input));
for i = 1:size(input, 1)
  for j = 1:size(input, 2)
    window = padded_input(i:i+kernel_size-1, j:j+kernel_size-1);
    output(i, j) = median(window(:));
  end
end
end