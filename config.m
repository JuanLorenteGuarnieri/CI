% name=['bunny_d=0.5_c=[256x256]',
%       'bunny_d=0.5_l=[1x1]_s=[256x256]',
%       'bunny_d=0.5_l=[16x16]_s=[16x16]',
%       'bunnybox_d=0.5_l=[16x16]_s=[16x16]',
%       'planes_d=0.5_l=[16x16]_s=[16x16]',
%       'Z_d=0.5_l=[1x1]_s=[256x256]'];
name = 'planes_d=0.5_l=[16x16]_s=[16x16]';             % Dataset name

%nVoxels=[2,4,8,16,32,64]
nVoxels = 64;                                     % Voxel resolution (x, y, z)

isFiltered = true;                              % Is the dataset filtered?
isAttenuated = true;                            % Is the dataset attenuated?
isPhasedFiltered = true;                        % Is the dataset phase filtered?

isPhasedOption1 = false;                         % Uses sigma as lambda_c/(2*log2), otherwise uses 2*lambda_c

dataset_path = strcat('data/', name, '.mat');           % Path to the dataset
volume_path = strcat('results/', name, '_', num2str(nVoxels), ...   % Path to save the volume
  '_', num2str(isFiltered), '_', num2str(isAttenuated), '_', num2str(isPhasedFiltered), '.mat');