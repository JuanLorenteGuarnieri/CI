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
filename = 'data/Z_d=0.5_l=[1x1]_s=[256x256].hdf5';  % <-- Update with your data file
dataset = load_hdf5_dataset(filename);

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

% Ensure volumePosition and volumeSize have 3 elements.
if numel(vol_center) < 3
  warning('Volume position has fewer than 3 elements. Filling missing dimensions with zeros.');
  vol_center = [vol_center; zeros(3-numel(vol_center),1)];
end
if numel(vol_size) < 3
  warning('Volume size has fewer than 3 elements. Repeating the provided value for missing dimensions.');
  vol_size = [vol_size; repmat(vol_size(1), 3-numel(vol_size), 1)];
end

voxel_res = [32, 32, 32];  % Voxel resolution (x, y, z)

% Create linearly spaced vectors for each dimension.
x_lin = linspace(vol_center(1) - vol_size(1)/2, vol_center(1) + vol_size(1)/2, voxel_res(1));
y_lin = linspace(vol_center(2) - vol_size(2)/2, vol_center(2) + vol_size(2)/2, voxel_res(2));
z_lin = linspace(vol_center(3) - vol_size(3)/2, vol_center(3) + vol_size(3)/2, voxel_res(3));

% Create the 3D grid of voxel centers.
[X, Y, Z] = ndgrid(x_lin, y_lin, z_lin);
num_voxels = numel(X);

% Initialize the reconstructed volume.
G = zeros(voxel_res);

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
  num_laser = size(laserPos,1);
  num_spad  = size(spadPos,1);
end

%% 4. Back-projection reconstruction
fprintf('Starting back-projection reconstruction...\n');
tic;
if isConfocal
  % Confocal data: iterate over each voxel and for each measurement position.
  for idx = 1:num_voxels
    xv = [X(idx), Y(idx), Z(idx)];  % current voxel center
    voxel_sum = 0;
    % Iterate over each measurement position (laser and SPAD co-located)
    for p = 1:num_positions
      pos = positions(p, :);
      % Calculate the 4 segments of the path:
      %   1. Laser device to the position on the wall
      %   2. Wall to voxel
      %   3. Voxel back to the same position on the wall
      %   4. Position on the wall to the SPAD device
      d1 = norm(origin_laser - pos);
      d2 = norm(pos - xv);
      d3 = norm(xv - pos);
      d4 = norm(pos - origin_spad);
      total_distance = d1 + d2 + d3 + d4;

      % Convert the total optical distance to a time index.
      t_bin = round((total_distance - t0)/deltaT) + 1;
      % Check that the index is within bounds.
      if t_bin >= 1 && t_bin <= size(H,3)
        voxel_sum = voxel_sum + H(p, t_bin);
      end
    end
    G(idx) = voxel_sum;
  end
else
  % Non-confocal data: iterate over each voxel and sum over laser-SPAD combinations.
  for idx = 1:num_voxels
    xv = [X(idx), Y(idx), Z(idx)];
    voxel_sum = 0;
    for i = 1:num_laser
      pos_laser = laserPos(i, :);
      for j = 1:num_spad
        pos_spad = spadPos(j, :);
        % Calculate the 4 segments of the path:
        %   1. From the laser device to the laser position on the wall
        %   2. From the laser position to the voxel
        %   3. From the voxel to the SPAD position on the wall
        %   4. From the SPAD position to the SPAD device
        d1 = norm(origin_laser - pos_laser);
        d2 = norm(pos_laser - xv);
        d3 = norm(xv - pos_spad);
        d4 = norm(pos_spad - origin_spad);
        total_distance = d1 + d2 + d3 + d4;

        t_bin = round((total_distance - t0)/deltaT) + 1;
        if t_bin >= 1 && t_bin <= size(H,3)
          voxel_sum = voxel_sum + H(i, j, t_bin);
        end
      end
    end
    G(idx) = voxel_sum;
  end
end
recon_time = toc;
fprintf('Reconstruction completed in %.2f seconds.\n', recon_time);

%% 5. Optional filtering: Apply Laplacian filter to enhance contrasts
f_lap = fspecial3('lap');
G_filtered = imfilter(G, -f_lap, 'symmetric');

%% 6. Visualization of the reconstructed volume
% Use volshow (available in recent MATLAB versions) with default settings.
figure;
volshow(G_filtered);
title('Reconstructed Volume with Laplacian Filtering');

% Alternatively, you can use volumeViewer:
% volumeViewer(G_filtered, 'Colormap', hot);

%% End of Script
%
% It is recommended to comment on parameter choices, reconstruction times, and
% make comparisons if testing with different datasets.
