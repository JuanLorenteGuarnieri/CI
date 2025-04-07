% run_experiments.m
% Script to execute lab4.m with different configurations

% Read the CSV file
data = readtable('results.csv', 'Delimiter', ',', 'ReadVariableNames', false);

% Assign column names for better readability
data.Properties.VariableNames = {'Name', 'nVoxels', 'isFiltered', ...
  'isAttenuated', 'isPhasedFiltered', 'isPhasedOption1', 'ExecutionTime'};

% Define the parameter values to test
names = {'bunny_d=0.5_c=[256x256]', 'bunny_d=0.5_l=[1x1]_s=[256x256]', ...
  'bunny_d=0.5_l=[16x16]_s=[16x16]', 'bunnybox_d=0.5_l=[16x16]_s=[16x16]', ...
  'planes_d=0.5_l=[16x16]_s=[16x16]'};
nVoxelsList = [2, 4, 8, 16, 32, 64];
isFilteredList = {true};
isAttenuatedList = {true};
isPhasedFilteredList = {true};
isPhasedOption1List = {true, false};
uniqueNames = unique(names);

for i = 1:length(nVoxelsList)
  nVoxels = nVoxelsList(i);
  for j = 1:length(uniqueNames)
    name = uniqueNames{j};
    for k = 1:length(isFilteredList)
      isFiltered = isFilteredList{k};
      for l = 1:length(isAttenuatedList)
        isAttenuated = isAttenuatedList{l};
        for m = 1:length(isPhasedFilteredList)
          isPhasedFiltered = isPhasedFilteredList{m};
          for n = 1:length(isPhasedOption1List)
            isPhasedOption1 = isPhasedOption1List{n};

            % Check if all parameters match a row in the data
            match = false;
            for row = 1:height(data)
              if strcmp(data.Name{row}, name) && ...
                  data.nVoxels(row) == nVoxels && ...
                  data.isFiltered(row) == double(isFiltered) && ...
                  data.isAttenuated(row) == double(isAttenuated) && ...
                  data.isPhasedFiltered(row) == double(isPhasedFiltered) && ...
                  data.isPhasedOption1(row) == double(isPhasedOption1)
                match = true;
                break;
              end
            end
            if ~match
              % Update the parameters in config.m without overwriting other content
              config_file = fileread('config.m');

              % Update or add the variables
              config_file = regexprep(config_file, 'name\s*=\s*''.*?'';', sprintf('name = ''%s'';', name));
              config_file = regexprep(config_file, 'nVoxels\s*=\s*\d+;', sprintf('nVoxels = %d;', nVoxels));
              config_file = regexprep(config_file, 'isFiltered\s*=\s*\w+;', sprintf('isFiltered = %s;', mat2str(isFiltered)));
              config_file = regexprep(config_file, 'isAttenuated\s*=\s*\w+;', sprintf('isAttenuated = %s;', mat2str(isAttenuated)));
              config_file = regexprep(config_file, 'isPhasedFiltered\s*=\s*\w+;', sprintf('isPhasedFiltered = %s;', mat2str(isPhasedFiltered)));
              config_file = regexprep(config_file, 'isPhasedOption1\s*=\s*\w+;', sprintf('isPhasedOption1 = %s;', mat2str(isPhasedOption1)));

              % Write the updated content back to config.m
              fid = fopen('config.m', 'w');
              fwrite(fid, config_file);
              fclose(fid);

              % Run the lab4.m script with dataset_path as input
              run('lab4.m');
            else
              fprintf('Parameters already exist in the data: %s, %d, %d, %d, %d, %d\n', ...
                name, nVoxels, isFiltered, isAttenuated, isPhasedFiltered, isPhasedOption1);
            end
          end
        end
      end
    end
  end
end