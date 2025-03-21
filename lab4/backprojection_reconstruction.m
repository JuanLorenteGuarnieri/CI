% backprojection_reconstruction.m
%
% Naive back-projection reconstruction for non-line-of-sight (NLOS) imaging.
% This script implements the algorithm described in the assignment and in Velten et al. (2012).
%
% Dependencies:
%   - load_hdf5_dataset.m (provided)
%   - fspecial3 (available on MATLAB File Exchange or implemented by the user)
%
% NOTE: Update the variable "filename" with the name/location of your dataset.
%
% Author: [Your Name]
% Date: [Today's Date]

clear; clc; close all;

%% 1. Load the dataset
filename = 'data/bunny_d=0.5_c=[256x256].mat';  % <-- Update with your data file
% filename = 'data/planes_d=0.5_l=[16x16]_s=[16x16].mat';
%dataset = load_hdf5_dataset(filename);
dataset = load(filename).data;

% Determine if the dataset is confocal
if ndims(dataset.data) == 3
  isConfocal = true;
  fprintf('The dataset is confocal.\n');
else
  isConfocal = false;
  fprintf('The dataset is non-confocal.\n');
end

%% 2. Define the voxel grid of the hidden scene
% Use volumePosition (center) and volumeSize (dimensions) from the dataset.
vol_center = dataset.volumePosition;   % 3x1 vector [x, y, z]
vol_size   = dataset.volumeSize;       % 3x1 vector [sx, sy, sz]

if isscalar(vol_size)
  vol_size = [vol_size, vol_size, vol_size];
end

num_voxels = [20,20,20];  % Voxel resolution (x, y, z)


% Compute the voxel spacing (resolution) in each direction
voxel_res = vol_size ./ num_voxels;
% Generate 1D coordinate grids for each axis
x = linspace(vol_center(1) - vol_size(1)/2 + voxel_res(1)/2, vol_center(1) + vol_size(1)/2 - voxel_res(1)/2, num_voxels(1));
y = linspace(vol_center(2) - vol_size(2)/2 + voxel_res(2)/2, vol_center(2) + vol_size(2)/2 - voxel_res(2)/2, num_voxels(2));
z = linspace(vol_center(3) - vol_size(3)/2 + voxel_res(3)/2, vol_center(3) + vol_size(3)/2 - voxel_res(3)/2, num_voxels(3));


% Initialize the reconstructed volume.
G = zeros(num_voxels);

%% 3. Extract measurement parameters and positions
t0 = dataset.t0;          % temporal offset (in meters, since deltaT is in optical distance units)
deltaT = dataset.deltaT;  % temporal resolution (in meters)
H = dataset.data;         % measurement data H (dimensions depend on confocal vs non-confocal)

if isConfocal
  % For confocal data, the laser and SPAD positions are the same.
  positions = dataset.spadPositions;  % Nx3 matrix
  origin_laser = dataset.laserOrigin; % 3x1 vector
  origin_spad  = dataset.spadOrigin;  % 3x1 vector
  num_positions = size(positions,1);
else
  % For non-confocal data, use separate grids for laser and SPAD.
  laserPos = dataset.laserPositions;  % Mx3 matrix
  spadPos  = dataset.spadPositions;   % Nx3 matrix
  origin_laser = dataset.laserOrigin;
  origin_spad  = dataset.spadOrigin;
end

%% 4. Back-projection reconstruction
fprintf('Starting back-projection reconstruction...\n');
tic;
if isConfocal
  % For confocal measurements, the laser and SPAD are at the same positions.
  for idx = 1:num_voxels(1)
    for idy = 1:num_voxels(2)
      for idz = 1:num_voxels(3)
        %xv = voxel_pos(idx, idy, idz);
        xv = [x(idx), y(idy), z(idz)];

        voxel_sum = 0;
        laser_size = size(positions);
        for i = 1:laser_size(1)
          for j = 1:laser_size(2)
            pos_laser = positions(i, j, :);
            pos_laser = reshape(pos_laser,1,[]); % or squeeze() transposed

            % 4 segments of the path:
            %   1. From the laser device to the laser position on the wall
            %   2. From the laser position to the voxel
            %   3. From the voxel back to the wall point (confocal: same as d2).
            %   4. From the SPAD position to the SPAD device
            d1 = norm(pos_laser(:) - origin_laser(:));
            d2 = norm(xv(:) - pos_laser(:));
            d3 = d2;%= norm(pos_laser(:) - xv(:));
            d4 = norm(origin_spad(:) - pos_laser(:));
            total_distance = d1 + d2 + d3 + d4;

            tv = round((total_distance - t0)/deltaT);
            if tv >= 1 && tv <= size(H,3)
              G(idx, idy, idz) = G(idx, idy, idz) + H(i, j, tv);
            end
          end
        end
      end
    end
  end
else
  % Non-confocal data: iterate over each voxel and sum over laser-SPAD combinations.
  for idx = 1:num_voxels(1)
    for idy = 1:num_voxels(2)
      for idz = 1:num_voxels(3)
        %xv = voxel_pos(idx, idy, idz);
        xv = [x(idx), y(idy), z(idz)];

        voxel_sum = 0;
        laser_size = size(laserPos);
        for i = 1:laser_size(1)
          for j = 1:laser_size(2)
            pos_laser = laserPos(i, j, :);
            pos_laser = reshape(pos_laser,1,[]); % or squeeze() transposed

            spadSize = size(spadPos);
            for xsi = 1:spadSize(1)
              for xsj = 1:spadSize(2)
                pos_spad = spadPos(xsi, xsj, :);
                pos_spad = reshape(pos_spad,1,[]); % or squeeze() transposed

                % 4 segments of the path:
                %   1. From the laser device to the laser position on the wall
                %   2. From the laser position to the voxel
                %   3. From the voxel to the SPAD position on the wall
                %   4. From the SPAD position to the SPAD device
                d1 = norm(pos_laser(:) - origin_laser(:));
                d2 = norm(xv(:) - pos_laser(:));
                d3 = norm(pos_spad(:) - xv(:));
                d4 = norm(origin_spad(:) - pos_spad(:));
                total_distance = d1 + d2 + d3 + d4;

                tv = round((total_distance - t0)/deltaT);
                if tv >= 1 && tv <= size(H,5)
                  G(idx, idy, idz) = G(idx, idy, idz) + H(i, j, xsi, xsj, tv);
                end
              end
            end
          end
        end
      end
    end
  end
end

recon_time = toc;
fprintf('Reconstruction completed in %.2f seconds.\n', recon_time);

%% 5. Optional filtering: Apply Laplacian filter to enhance contrasts
f_lap = fspecial3('lap');
%G_filtered = G;
G_filtered = imfilter(G, -f_lap, 'symmetric');

%% 6. Visualization of the reconstructed volume
% Use volshow (available in recent MATLAB versions) with default settings.
figure;
vs_h = volshow(G_filtered, RenderingStyle="MaximumIntensityProjection", Colormap=hot);
vs_h.Parent.BackgroundColor = [0 0 0];
vs_h.Parent.GradientColor = [0 0 0];
title('Reconstructed Volume with Laplacian Filtering');

% Alternatively, you can use volumeViewer:
%volumeViewer(G_filtered, 'Colormap', hot);

%% End of Script
%
% It is recommended to comment on parameter choices, reconstruction times, and
% make comparisons if testing with different datasets.
