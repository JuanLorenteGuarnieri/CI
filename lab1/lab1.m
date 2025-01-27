%%%%%%%%%%%%%%%%%%%     2.1 Reading the image into Matlab     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the image in TIFF format
workspacePath = fileparts(mfilename('fullpath'));
imgPath = fullfile(workspacePath, '/images_tiff/IMG_1026.tiff');
img = imread(imgPath);

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
figure;
imshow(img);
title('Original Image');

% Display the double-precision image
figure;
imshow(img_double, []);
title('Double-Precision Image');

%%%%%%%%%%%%%%%%%%%             2.2 Linearization             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%               3 Demosaicing               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%             4 White balancing             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%                5 Denoising                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%              6 Color balance              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%            7 Tone reproduction            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%               8 Compression               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


