% Config variables
color = true;
method = 'zDeconvWNR'; % 'zDeconvWNR' or 'deconvwnr' or 'deconvlucy'
aperture_method = 'circular'; % 'circular' or 'Levin' or 'raskar' or 'zhou' or 'rand'
percentage_rand_aperture = 0.9;

sigma = 0.5; % Noise level (Gaussian noise)
blurSize = 7; % Blur size
iter = 150; % Number of iterations for deconvlucy


% Read data
workspacePath = fileparts(mfilename('fullpath'));
if strcmp(aperture_method, 'rand')
    matrix = rand(100,100) < percentage_rand_aperture;
    aperture = uint8(matrix* 255);
else
    aperture = imread(fullfile(workspacePath, '/apertures/', [aperture_method, '.bmp']));
end
image = imread(fullfile(workspacePath, '/images/', 'burano.jpg'));

if color == false
    image = image(:, :, 1);
end

f0 = im2double(image);
[height, width, channel] = size(f0);

% Prior matrix: 1/f law
A_star = eMakePrior(height, width) + 0.00000001;
C = sigma.^2 * height * width ./ A_star;

% Normalization
temp = fspecial('disk', blurSize);
flow = max(temp(:));

% Calculate effective PSF
k1 = im2double(...
    imresize(aperture, [2*blurSize + 1, 2*blurSize + 1], 'nearest')...
    );

k1 = k1 * (flow / max(k1(:)));

if color
    % Apply blur and recover each channel independently
    f1 = zeros(size(f0));
    f0_hat = zeros(size(f0));

    for c = 1:channel
        % Apply blur
        f1(:, :, c) = zDefocused(f0(:, :, c), k1, sigma, 0);

        % Recover
        if strcmp(method, "zDeconvWNR")
            f0_hat(:, :, c) = zDeconvWNR(f1(:, :, c), k1, C);
        end
        if strcmp(method, "deconvwnr")
            f0_hat(:, :, c) = deconvwnr(f1(:, :, c), k1, C); % deconvwnr(I,psf,nsr)
        end
        if strcmp(method, "deconvlucy")
            f0_hat(:, :, c) = deconvlucy(f1(:, :, c), k1, iter); % deconvlucy(I,psf,iter)
        end
    end
else
    % Apply blur
    f1 = zDefocused(f0, k1, sigma, 0);

    % Recover
    if strcmp(method, "zDeconvWNR")
        f0_hat = zDeconvWNR(f1, k1, C);
    end
    if strcmp(method, "deconvwnr")
        f0_hat = deconvwnr(f1, k1, C); % deconvwnr(I,psf,nsr)
    end
    if strcmp(method, "deconvlucy")
        f0_hat = deconvlucy(f1, k1, iter); % deconvlucy(I,psf,iter)
    end

end


% Display results
figure;

subplot_tight(1, 3, 1, 0.0, false)
imshow(f0);
title('Focused');

subplot_tight(1, 3, 2, 0.0, false)
imshow(f1);
title('Defocused');

subplot_tight(1, 3, 3, 0.0, false)
imshow(f0_hat);
title('Recovered');
