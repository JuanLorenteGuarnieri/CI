% Script to read and process the CSV file

% Define the file path
filePath = 'results.csv';

% Define the execution type
execution_type = 'part3'; % Options: 'individual', 'group', 'part2', 'part3'

% Read the CSV file
data = readtable(filePath, 'Delimiter', ',', 'ReadVariableNames', false);

% Assign column names for better readability
data.Properties.VariableNames = {'Name', 'nVoxels', 'isFiltered', ...
  'isAttenuated', 'isPhasedFiltered', 'isPhasedOption1', 'ExecutionTime'};

% Get unique names
uniqueNames = unique(data.Name);
if strcmp(execution_type, 'individual')
  % Loop through each unique name
  for i = 1:length(uniqueNames)
    % Filter data for the current name
    groupData = data(strcmp(data.Name, uniqueNames{i}), :);

    % Create a new figure for the current name
    figure;
    hold on;

    % Define a colormap for the subgroups
    colors = lines(3); % Assuming 3 binary variables

    % Plot data for each combination of the other variables
    for isFiltered = 0:1
      for isAttenuated = 0:1
        for isPhasedFiltered = 0:1
          % Filter rows for the current combination of variables
          subgroupData = groupData(groupData.isFiltered == isFiltered & ...
            groupData.isAttenuated == isAttenuated & ...
            groupData.isPhasedFiltered == isPhasedFiltered, :);

          % Skip if no data for this combination
          if isempty(subgroupData)
            continue;
          end

          % Plot the data
          plot(subgroupData.nVoxels, subgroupData.ExecutionTime, 'o-', ...
            'Color', colors(isFiltered + isAttenuated + 1, :), ...
            'DisplayName', sprintf('Filtered=%d, Attenuated=%d, Phased=%d', ...
            isFiltered, isAttenuated, isPhasedFiltered));
        end
      end
    end

    % Set logarithmic scale for the y-axis
    set(gca, 'YScale', 'log');

    % Define custom ticks for the y-axis
    customTicks = [0, 1, 5, 10, 50, 100, 500, 1000, 5000]; % Adjust values as needed
    set(gca, 'YTick', customTicks);

    % Add labels, title, grid, and legend
    xlabel('Number of Voxels');
    ylabel('Execution Time (s)');
    title(['Execution Time for Name: ', uniqueNames{i}]);
    grid on;
    legend('show', 'Location', 'southeast');
    hold off;

  end
elseif strcmp(execution_type, 'part3')
  % Loop through each unique name starting with 'bunny'
  for i = 1:length(uniqueNames)
    if startsWith(uniqueNames{i}, 'bunny')
      % Filter data for the current name
      groupData = data(strcmp(data.Name, uniqueNames{i}), :);

      % Create a new figure for the current name
      figure;
      hold on;

      % Define a colormap for the subgroups
      colors = lines(2); % Assuming only 2 states for isFiltered (0 or 1)

      % Plot data for each filtering state
      for isFiltered = 0:1
        % Filter rows for the current filtering state
        subgroupData = groupData(groupData.isFiltered == isFiltered, :);

        % Skip if no data for this filtering state
        if isempty(subgroupData)
          continue;
        end

        % Plot the data
        plot(subgroupData.nVoxels, subgroupData.ExecutionTime, 'o-', ...
          'Color', colors(isFiltered + 1, :), ...
          'DisplayName', sprintf('Filtered=%d', isFiltered));
      end

      % Set logarithmic scale for the y-axis
      set(gca, 'YScale', 'log');

      % Define custom ticks for the y-axis
      customTicks = [0, 1, 5, 10, 50, 100, 500, 1000, 5000]; % Adjust values as needed
      set(gca, 'YTick', customTicks);

      % Add labels, title, grid, and legend
      xlabel('Number of Voxels');
      ylabel('Execution Time (s)');
      title(['Execution Time for Name: ', uniqueNames{i}]);
      grid on;
      legend('show', 'Location', 'southeast');
      hold off;
    end
  end
elseif strcmp(execution_type, 'group')
  % Plot all data for the current name as a single group
  plot(groupData.nVoxels, groupData.ExecutionTime, 'o-', ...
    'DisplayName', 'Group Data');
elseif strcmp(execution_type, 'part2')
  % Filter for names starting with 'Z' and isFiltered == true
  for i = 1:length(uniqueNames)
    % Filter data for the current name
    if startsWith(uniqueNames{i}, 'Z')
      subgroupData = data(strcmp(data.Name, uniqueNames{i}) & data.isFiltered == true & data.isAttenuated == false, :);

      % Create a new figure for the current name
      figure;
      hold on;

      % Define a colormap for the subgroups
      colors = lines(3); % Assuming 3 binary variables


      % Skip if no data for this combination
      if isempty(subgroupData)
        continue;
      end

      % Plot the data
      % plot(subgroupData.nVoxels, subgroupData.ExecutionTime, 'o-', ...
      %   'Color', colors(isFiltered + isAttenuated + 1, :), ...
      %   'DisplayName', sprintf('Filtered=%d, Attenuated=%d, Phased=%d', ...
      %   isFiltered, isAttenuated, isPhasedFiltered));
      % Plot the data
      plot(subgroupData.nVoxels, subgroupData.ExecutionTime, 'o-', ...
        'DisplayName', ['Filtered Data for ', uniqueNames{i}]);

      % Set logarithmic scale for the y-axis
      set(gca, 'YScale', 'log');

      % Define custom ticks for the y-axis
      customTicks = [0, 1, 5, 10, 50, 100, 500, 1000, 5000]; % Adjust values as needed
      set(gca, 'YTick', customTicks);

      % Add labels, title, grid, and legend
      xlabel('Number of Voxels');
      ylabel('Execution Time (s)');
      title(['Execution Time for Name: ', uniqueNames{i}]);
      grid on;
      legend('show', 'Location', 'southeast');
      hold off;
    end
  end
else
  error('Invalid execution_type. Use "individual", "group", or "part2".');
end