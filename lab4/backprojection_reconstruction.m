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
filename = 'data/Z_d=0.5_l=[1x1]_s=[256x256].mat';  % <-- Update with your data file
%filename = 'data/planes_d=0.5_l=[16x16]_s=[16x16].mat';
%dataset = load_hdf5_dataset(filename);
dataset = load(filename).data;

% Determine if the dataset is confocal
if ndims(dataset.data) == 3
  isConfocal = t13rue;
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

% Ensure volumePosition and volumeSize have 3 elements.
%if numel(vol_center) < 3
%  warning('Volume position has fewer than 3 elements. Filling missing dimensions with zeros.');
%  vol_center = [vol_center; zeros(3-numel(vol_center),1)];
%end
%if numel(vol_size) < 3
%  warning('Volume size has fewer than 3 elements. Repeating the provided value for missing dimensions.');
%  vol_size = [vol_size; repmat(vol_size(1), 3-numel(vol_size), 1)];
%end

voxel_res = [8, 8, 8];  % Voxel resolution (x, y, z)

% Create linearly spaced vectors for each dimension.
x_lin = linspace(vol_center(1) - vol_size(1)/2, vol_center(1) + vol_size(1)/2, voxel_res(1));
y_lin = linspace(vol_center(2) - vol_size(2)/2, vol_center(2) + vol_size(2)/2, voxel_res(2));
z_lin = linspace(vol_center(3) - vol_size(3)/2, vol_center(3) + vol_size(3)/2, voxel_res(3));

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
end

%% 4. Back-projection reconstruction
fprintf('Starting back-projection reconstruction...\n');
tic;
if isConfocal
  % Confocal data: iterate over each voxel and for each measurement position.
  %for idx = 1:num_voxels
  %  xv = [X(idx), Y(idx), Z(idx)];  % current voxel center
  %  voxel_sum = 0;
  %  % Iterate over each measurement position (laser and SPAD co-located)
  %  for p = 1:num_positions
  %    pos = positions(p, :);
  %    % Calculate the 4 segments of the path:
  %    %   1. Laser device to the position on the wall
  %    %   2. Wall to voxel
  %    %   3. Voxel back to the same position on the wall
  %    %   4. Position on the wall to the SPAD device
  %    d1 = norm(pos - origin_laser);
  %    d2 = norm(xv - pos);
  %    d3 = norm(pos - xv);
  %    d4 = norm(origin_spad - pos);
  %    total_distance = d1 + d2 + d3 + d4;

  %    % Convert the total optical distance to a time index.
  %    t_bin = round((total_distance - t0)/deltaT) + 1;
  %    % Check that the index is within bounds.
  %    if t_bin >= 1 && t_bin <= size(H,3)
  %      voxel_sum = voxel_sum + H(p, t_bin);
  %    end
  %  end
  %  G(idx) = voxel_sum;
  %end
else
  % Non-confocal data: iterate over each voxel and sum over laser-SPAD combinations.
  for idx = 1:length(x_lin)
  for idy = 1:length(y_lin)
  for idz = 1:length(z_lin)
    xv = [x_lin(idx), y_lin(idy), z_lin(idz)];
    voxel_sum = 0;
    laser_size = size(laserPos);
    for i = 1:laser_size(1)
    for j = 1:laser_size(2)
      pos_laser = laserPos(i, j, :);
      pos_laser = reshape(pos_laser,1,[]); % or squeeze() transposed

      spadSize = size(spadPos);
      for xsi = 1:4:spadSize(1)
      for xsj = 1:4:spadSize(2)
        pos_spad = spadPos(xsi, xsj, :);
        pos_spad = reshape(pos_spad,1,[]); % or squeeze() transposed

        % 4 segments of the path:
        %   1. From the laser device to the laser position on the wall
        %   2. From the laser position to the voxel
        %   3. From the voxel to the SPAD position on the wall
        %   4. From the SPAD position to the SPAD device
        d1 = norm(pos_laser - origin_laser);
        d2 = norm(xv - pos_laser);
        d3 = norm(pos_spad - xv);
        d4 = norm(origin_spad - pos_spad);
        total_distance = d1 + d2 + d3 + d4;
        
        t_bin = round((total_distance - t0)/deltaT) + 1;
        if t_bin >= 1 && t_bin <= size(H,5)
          voxel_sum = voxel_sum + H(i, j, xsi, xsj, t_bin);
        end
      end
      end
    end
    end
    
    G(idx, idy, idz) = voxel_sum; 
  end  
  end
  end
end

recon_time = toc;
fprintf('Reconstruction completed in %.2f seconds.\n', recon_time);

%% 5. Optional filtering: Apply Laplacian filter to enhance contrasts
f_lap = fspecial3('lap');
G_filtered = G;
%G_filtered = imfilter(G, -f_lap, 'symmetric');

%% 6. Visualization of the reconstructed volume
% Use volshow (available in recent MATLAB versions) with default settings.
figure;
vs_h = volshow(G, RenderingStyle="MaximumIntensityProjection", Colormap=hot);
vs_h.Parent.BackgroundColor = [0 0 0];
vs_h.Parent.GradientColor = [0 0 0];
title('Reconstructed Volume with Laplacian Filtering');

% Alternatively, you can use volumeViewer:
%volumeViewer(G_filtered, 'Colormap', hot);

%% End of Script
%
% It is recommended to comment on parameter choices, reconstruction times, and
% make comparisons if testing with different datasets.
