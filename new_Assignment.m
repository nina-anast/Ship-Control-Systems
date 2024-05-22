%% *   *
% 
%% * 1:*         
% 

clear all
warning("off")
% load all data
Data2024 = readtable("Data2024.txt",'TrailingDelimitersRule','ignore')
% Extract the units
% Open the file for reading
fid = fopen('Data2024.txt', 'r');
% Read the header line and discard it
fgetl(fid);
% Read the next line as a string
second_line = fgetl(fid);
% Close the file
fclose(fid);
% Split the string into cell array of units
units = strsplit(second_line, '\t');
% Convert cell array to array
units = string(units);

clearvars fid second_line ans
%%
% Initialize the fieldnames of the data
fieldnames = Data2024.Properties.VariableNames
% Initialize a struct to store values
data = struct();

% Loop through each column
for i = 1:numel(fieldnames)
    % Extract values of the current column
    data.(fieldnames{i}) = table2array(Data2024(:, fieldnames{i}));
end

% Define the number of plots
num_plots = numel(fieldnames) - 1; % Excluding the first column

% Define the number of rows and columns for subplots
num_rows = ceil(sqrt(num_plots));
num_cols = ceil(num_plots / num_rows);

% Create a new figure
figure;

% Plot each column against the "Time" column on subplots
for i = 2:numel(fieldnames) % Start from the second column
    subplot(num_rows, num_cols, i-1); % Create subplot
    plot(data.Time, data.(fieldnames{i}),'.');
    title(fieldnames{i}); % Set title as column name
    xlabel('Time [s]'); % X-axis label
    ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Y-axis label
end
hold off
clearvars num_cols num_rows num_plots i ans
%  

% i = 2;
% figure; % Create a new figure
% plot(data.Time, data.(fieldnames{i}),".");
% title(fieldnames{i}); % Set title as column name
% xlabel('Time [s]'); % X-axis label
% 
% % Concatenate field name with unit and set it as Y-axis label
% ylabel(sprintf('%s %s', fieldnames{i}, units{i}));
% clearvars i
%% 
% *NOx:*     , .   225s. 
%         ,   600 ppm   
%   200 ppm.
% 
% *Fuel Consumption:*       (  
%  x 180-405s   ).      20kg/h  
%   65kg/h.
% 
% *:*     20    2619s    
%  peak.   8 peaks        130.
%%  2:     .
%           , 
%      threshold.

% % Initialize a struct to store filtered values
% data_filt = struct();
% time_filt = struct();
% 
% % Store thresholds
% thresholds = zeros(1,length(fieldnames));
% 
% clc
% 
% for j = 2:numel(fieldnames)
%     [filteredData, filteredTime,threshold] = chauvenetFilter(Data2024,fieldnames,units,j);
%     data_filt.(fieldnames{j}) = filteredData;
%     time_filt.(fieldnames{j}) = filteredTime;
%     thresholds(1,j) = threshold;
% end
% clearvars j threshold filteredData filteredTime
%   
%          thresholds  
%          .

% writematrix(thresholds,"thresholds.csv");
% clearvars ans 
%  
%      ,     
%  :

% % Initialize a struct to store filtered values
% data_filt = struct();
% time_filt = struct();
% 
% % Read thresholds from the CSV file
% thresholds = readmatrix('thresholds.csv');
% 
% for j = 2:numel(fieldnames)
%     [filteredData, filteredTime, ~] = chauvenetFilter2(Data2024, fieldnames, units, j, thresholds(j));
%     data_filt.(fieldnames{j}) = filteredData;
%     time_filt.(fieldnames{j}) = filteredTime;
% end
% 
% clearvars filteredData filteredTime j
%% 
%  *Interquartile Range*:
% 
% Apply to:
%% 
% * Fuel Consumption (3)
% * Intake Pressure (6)
% * Torque Reference (7)
% * Rot Speed (8)
% * Engine Torque (9)
% * (EGR Command(10) ->  threshold)
% * Exhaust Gas Temperature (11)

% % Initialize a struct to store filtered values
% data_filt = struct();
% time_filt = struct();
% 
% for j = 2:numel(fieldnames)
%     [filteredData,filteredTime] = IQR(Data2024.(fieldnames{j}),Data2024.Time);
%     data_filt.(fieldnames{j}) = filteredData;
%     time_filt.(fieldnames{j}) = filteredTime;
% end 
% 
% clearvars filteredData filteredTime j
%% 
% Plot for *filtered data*

% i = 3;
% figure; % Create a new figure
% plot(time_filt.(fieldnames{i}),data_filt.(fieldnames{i}),".");
% title(fieldnames{i}); % Set title as column name
% xlabel('Time [s]'); % X-axis label
% 
% % Concatenate field name with unit and set it as Y-axis label
% ylabel(sprintf('%s %s', fieldnames{i}, units{i}));
%% 
% *Remove outliers using median:*
% 
% use in:
%% 
% * fuel consumption (3)
% * exhaust mass flow (5)
% * Torque Reference (7)
% * Rot Speed (8)
% * Engine Torque (9)
% * Exhaust Gas Temperature (11)

% % Initialize a struct to store filtered values and corresponding times
% data_filt = struct();
% time_filt = struct();
% 
% i = 11;
% % Filter the data using rmoutliers
% [data_filt.(fieldnames{i}), idx] = rmoutliers(Data2024.(fieldnames{i}), "median");
% 
% % Retain corresponding times using logical indexing
% time_filt.(fieldnames{i}) = Data2024.Time;
% time_filt.(fieldnames{i})(idx) = [];
% 
% % Plot unfiltered data
% figure;
% subplot(2,1,1); % Plot on the first row, first two columns
% plot(Data2024.Time, Data2024.(fieldnames{i}),".");
% title(fieldnames{i}); % Set title as column name
% xlabel('Time [s]'); % X-axis label
% ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label
% 
% % Plot unfiltered data
% subplot(2,1,2); % Plot on the second row, first column
% plot(time_filt.(fieldnames{i}),data_filt.(fieldnames{i}),".");
% title(fieldnames{i}); % Set title as column name
% xlabel('Time [s]'); % X-axis label
% ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label
%% 
% * :*

% % Initialize a struct to store filtered values
% data_filt = struct();
% time_filt = struct();
% 
% % Read thresholds from the CSV file
% thresholds = readmatrix('thresholds.csv');
% 
% for j = 2:numel(fieldnames)
%     if j == 2 || j == 4 || j == 10
%         [filteredData, filteredTime] = chauvenetFilter2(Data2024, fieldnames, units, j, thresholds(j));
%         % Calculate removedIndices based on original and filtered data
%         removedIndices = find(~ismember(Data2024.Time, filteredTime));
%         data_filt.(fieldnames{j}) = filteredData;
%         time_filt.(fieldnames{j}) = filteredTime;
%         removedTimes = unique(Data2024.Time(removedIndices));
%         time_filt.removed.(fieldnames{j}) = removedTimes;
%         % Remove corresponding data
%         Data2024(removedIndices,:) = [];
%     elseif j == 3 || j == 6 || j == 7
%         [filteredData, filteredTime] = IQR(Data2024.(fieldnames{j}), Data2024.Time);
%         % Calculate removedIndices based on original and filtered data
%         removedIndices = find(~ismember(Data2024.Time, filteredTime));
%         data_filt.(fieldnames{j}) = filteredData;
%         time_filt.(fieldnames{j}) = filteredTime;
%         removedTimes = unique(Data2024.Time(removedIndices));
%         time_filt.removed.(fieldnames{j}) = removedTimes;
%         % Remove corresponding data
%         Data2024(removedIndices,:) = [];
%     elseif j == 5 || j == 8 || j == 9 || j == 11
%         % Filter the data using rmoutliers
%         [data_filt.(fieldnames{j}), idx] = rmoutliers(Data2024.(fieldnames{j}), "median");
%         % Retain corresponding times using logical indexing
%         time_filt.(fieldnames{j}) = Data2024.Time;
%         removedIndices = idx;
%         removedTimes = unique(Data2024.Time(removedIndices));
%         time_filt.removed.(fieldnames{j}) = removedTimes;
%         % Remove corresponding data
%         Data2024(removedIndices,:) = [];
%     end
% end
% 
% clearvars data_filt filteredData filteredTime idx j removedIndices removedTimes thresholds time_filt
%%
% % Save training and testing structs to files
% filtered_data_file = 'filtered_data.mat';
% save(filtered_data_file, 'Data2024');
% % plot_differences(Data2024,fieldnames,units,data_filt,time_filt,4)
% * :*
% 

% i = 8;
% fieldname = fieldnames{i};
% % Normalize Data
% [normalizedValidating,centerValueValidating,scaleValueValidating] = normalize(validating,"range");
% 
% % Display results
% figure
% tiledlayout(2,1);
% nexttile
% plot(filtered_data.(fieldnames{i}),"Color",[77 190 238]/255,"DisplayName","Input data")
% legend
% ylabel(fieldnames{i})
% 
% nexttile
% plot(normalizedValidating.(fieldnames{i}),"Color",[0 114 189]/255,"LineWidth",1.5,...
%     "DisplayName","Normalized data")
% legend
% ylabel(fieldnames{i})
% set(gcf,"NextPlot","New")
%  :

% % Initialize structs to store training and testing data
% training = struct();
% testing = struct();
% 
% for i = 1:numel(fieldnames)
%     % Define the percentage of data for training and testing
%     train_percent = 0.8;
%     test_percent = 1 - train_percent;
% 
%     % Get the field name
%     fieldname = fieldnames{i};
% 
%     % Extract the data directly from Data2024
%     data2 = normalizedTable.(fieldname);
%     times2 = normalizedTable.Time;
% 
%     % Create indices for cross-validation with 80% training and 20% testing
%     cv = cvpartition(size(data2, 1), 'Holdout', test_percent);
% 
%     % Get indices for training and testing sets
%     train_idx = cv.training;
%     test_idx = cv.test;
% 
%     % Split the data into training and testing sets
%     training.(fieldname) = data2(train_idx, :);
%     testing.(fieldname) = data2(test_idx, :);
% end
% 
% % Define file paths for saving
% training_file = 'training_data-normalized.mat';  % Specify the file path for saving training data
% testing_file = 'testing_data-normalized.mat';    % Specify the file path for saving testing data
% 
% % Save training and testing structs to files
% save(training_file, 'training');
% save(testing_file, 'testing');
%%
% % Initialize structs to store training and validating data
% training2 = struct();
% validating = struct();
% 
% for i = 1:numel(fieldnames)
%     % Define the percentage of data for training and testing
%     train_percent = 0.8;
%     val_percent = 1 - train_percent;
% 
%     % Get the field name
%     fieldname = fieldnames{i};
% 
%     % Extract the data and corresponding times corresponding to the field name
%     data2 = training.(fieldname);
%     times2 = training.Time;
% 
%     % Create indices for cross-validation with 80% training and 20% testing
%     cv = cvpartition(size(data2, 1), 'Holdout', val_percent);
% 
%     % Get indices for training and testing sets
%     train_idx = cv.training;
%     val_idx = cv.test;
% 
%     % Split the data into training and testing sets
%     training2.(fieldname) = data2(train_idx, :);
%     validating.(fieldname) = data2(val_idx, :);
% end
% 
% % Define file paths for saving
% training_file = 'training2_data-normalized.mat';  % Specify the file path for saving training data
% validating_file = 'validating_data-normalized.mat';    % Specify the file path for saving testing data
% 
% % Save training and validating structs to files
% save(training_file, 'training2');
% save(validating_file, 'validating');

%% 
% *To check the separation was ok:*

% i = 5;
% 
% % Plot unfiltered data
% figure;
% subplot(2,1,1); % Plot on the first row, first two columns
% plot(Data2024.Time, Data2024.(fieldnames{i}), ".");
% title(fieldnames{i}); % Set title as column name
% xlabel('Time [s]'); % X-axis label
% ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label
% 
% % Plot filtered data
% subplot(2,1,2); % Plot on the second row, first column
% plot(testing.Time, testing.(fieldnames{i}), ".");
% title(['Filtered ' fieldnames{i}]);
% xlabel('Time [s]'); % X-axis label
% ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label
% clearvars i
%% 
% *Load previewsly made training and testing data:*

% % Load previewsly filtered Data
% filtered_data = load("separated_data\filtered\filtered_data.mat");
% filtered_data = filtered_data.Data2024
% % Load the training data from the .mat file
% training = load('separated_data\filtered\training2_data.mat');
% training = struct2table(training.training2)
% % Load the testing data from the .mat file
% testing = load('separated_data\filtered\testing_data.mat');
% testing = struct2table(testing.testing)
% % Load the validating data from the .mat file
% validating = load('separated_data\filtered\validating_data.mat');
% validating = struct2table(validating.validating)
%% 
% *Load previewsly made training and testing data:*

% Load previewsly filtered Data
filtered_data = load("normalized_data.mat");
filtered_data = filtered_data.normalizedTable
% filtered_data = filtered_data.Data2024
% Load the training data from the .mat file
training = load('training_data-normalized.mat');
training = training.normalizedTraining
% Load the testing data from the .mat file
testing = load('testing_data-normalized.mat');
testing = testing.normalizedTesting
% Load the validating data from the .mat file
validating = load('validating_data-normalized.mat');
validating = validating.normalizedValidating
% Correlation Matrix

% % Convert the table to a matrix
% dataMatrix = table2array(filtered_data);
% correlationMatrix = corr(dataMatrix);
% figure;
% 
% % Plot correlation matrix
% % corrplot(dataMatrix, 'varNames', filtered_data.Properties.VariableNames);
%  
%         matrices

% Convert Time, RotSpeed and lambda's test, train and validation sets from tables to matrices
input1_train = table2array([training(:,1) training(:,4) training(:,8)]);
input1_val = table2array([validating(:,1) validating(:,4) validating(:,8)]);
input1_test = table2array([testing(:,1) testing(:,4) testing(:,8)]);

% Convert NOx's test, train and validation sets from tables to matrices
NOx_train = table2array(training(:,2));
NOx_val =  table2array(validating(:,2));
NOx_test = table2array(testing(:,2));

% Convert Time,lambda,ExhaustGasMassflow and Intake Prssure's test, train and validation sets from tables to matrices
input2_train = table2array([training(:,1) training(:,4) training(:,5) training(:,6)]);
input2_val = table2array([validating(:,1) validating(:,4) validating(:,5) validating(:,6)]);
input2_test = table2array([testing(:,1) testing(:,4) testing(:,5) testing(:,6)]);

% Convert Fuel Consumption's test, train and validation sets from tables to matrices
Fuel_Consumption_train = table2array(training(:,2));
Fuel_Consumption_val =  table2array(validating(:,2));
Fuel_Consumption_test = table2array(testing(:,2));
%%
% Input NOx network's archtecture as recommended
layers1 = [
    featureInputLayer(3)
    fullyConnectedLayer(100)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];

% Input NOx network's archtecture as recommended
layers2 = [
    featureInputLayer(4)
    fullyConnectedLayer(100)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];

% Open DeepNetworkDesigner App
deepNetworkDesigner

% Define training options 
options1 = trainingOptions('adam',...
    'Shuffle','every-epoch',...
    'MiniBatchSize',1000, ...
    'MaxEpochs',10, ...
    'InitialLearnRate',1e-3, ...
    'ValidationData',{input1_val, NOx_val},...
    'Plots','training-progress', ...
    'Verbose',false);


% Train 1st neural network to predict NOx emissions
NOx_net = trainNetwork(input1_train,NOx_train,layers1,options1);
%%
% Plot NOx emissions prediction using lambda's training, validation and test sets 

pred_and_plot(input1_train,NOx_train,NOx_net,1,2,fieldnames)
pred_and_plot(input1_val,NOx_val,NOx_net,1,2,fieldnames)
pred_and_plot2(input1_test,NOx_test,NOx_net,1,2,fieldnames,units,centerValueTesting,scaleValueTesting)
%%
% Define training options 
options2 = trainingOptions('adam',...
    'Shuffle','every-epoch',...
    'MiniBatchSize',1000, ...
    'MaxEpochs',10, ...
    'InitialLearnRate',1e-3, ...
    'ValidationData',{input2_val, Fuel_Consumption_val},...
    'Plots','training-progress', ...
    'Verbose',false);

% Train 2nd neural network to predict Fuel Consumption
Fuel_Consumption_net = trainNetwork(input2_train,Fuel_Consumption_train,layers2,options2);
%%
% Plot NOx emissions prediction using lambda's training, validation and test sets 
pred_and_plot(input2_train,Fuel_Consumption_train,Fuel_Consumption_net,1,3,fieldnames)
pred_and_plot(input2_val,Fuel_Consumption_val,Fuel_Consumption_net,1,3,fieldnames)
pred_and_plot2(input2_test,Fuel_Consumption_test,Fuel_Consumption_net,1,3,fieldnames,units,centerValueTesting,scaleValueTesting)
%%
% % *Get initial data from normalized:*
% 
% % Reverting the normalization
% originalTesting = testing .* scaleValueTesting + centerValueTesting;
% 
% % Display results
% figure
% tiledlayout(2,1);
% nexttile
% plot(originalTesting.(fieldnames{i}),"Color",[77 190 238]/255,"DisplayName","Input data")
% legend
% ylabel(fieldnames{i})
% 
% nexttile
% plot(normalizedTesting.(fieldnames{i}),"Color",[0 114 189]/255,"LineWidth",1.5,...
%     "DisplayName","Normalized data")
% legend
% ylabel(fieldnames{i})
% set(gcf,"NextPlot","New")
% Functions
%     *thresholds*:

function [filteredData, filteredTime, threshold] = chauvenetFilter(Data, fieldnames, units, i)

    data = Data.(fieldnames{i});
    time = Data.Time;
    
    redo = true;
    while redo
        % Calculate mean and standard deviation
        mu = mean(data);
        sigma = std(data);

        % Calculate residuals
        residuals = abs(data - mu);

        % Calculate probability
        p = 1 - 0.5 * (1 + erf(residuals / (sigma * sqrt(2))));
        p_mean = mean(p);
        
        % Plot unfiltered data
        figure;
        subplot(2,2,1); % Plot on the first row, first two columns
        plot(time, data);
        title(fieldnames{i}); % Set title as column name
        xlabel('Time [s]'); % X-axis label
        ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label

        % Plot unfiltered data
        subplot(2,2,2); % Plot on the second row, first column
        plot(time, data,".");
        title(fieldnames{i}); % Set title as column name
        xlabel('Time [s]'); % X-axis label
        ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label

        % Prompt user for threshold
        threshold = input(['Enter the threshold probability (mean is ', num2str(p_mean), '): ']);

        % Apply Chauvenet's Criterion (threshold probability)
        outliers = p < threshold;

        % Remove outliers
        filteredData = data(~outliers);
        filteredTime = time(~outliers);

        % Plot filtered data
        subplot(2,2,[3,4]); % Plot on the second row, second column
        plot(filteredTime, filteredData);
        title(['Filtered ', fieldnames{i}]); % Set title as column name
        xlabel('Time [s]'); % X-axis label
        ylabel(sprintf('%s %s (Filtered)', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label
        
        % Ask if the user wants to redo
        redo_prompt = input('Do you want to redo the filtering? (yes/no): ', 's');
        if strcmpi(redo_prompt, 'no')
            redo = false;
        end
    end
end
%% 
%   **   :

function [filteredData, filteredTime, threshold] = chauvenetFilter2(Data, fieldnames, units, i, threshold)

    data = Data.(fieldnames{i});
    time = Data.Time;
    
    % Calculate mean and standard deviation
    mu = mean(data);
    sigma = std(data);

    % Calculate residuals
    residuals = abs(data - mu);

    % Calculate probability
    p = 1 - 0.5 * (1 + erf(residuals / (sigma * sqrt(2))));
    p_mean = mean(p);

    % Apply Chauvenet's Criterion (threshold probability)
    outliers = p < threshold;

    % Remove outliers
    filteredData = data(~outliers);
    filteredTime = time(~outliers);
end
%% 
%  *Interquartile Range*:
% 
% <https://www.indeed.com/career-advice/career-development/outliers-statistics 
% https://www.indeed.com/career-advice/career-development/outliers-statistics>

function [FilteredData,Time]=IQR(Data,Time)
%% Remove outliers from vectors with interquartile range method

% Sort Input Data in ascending order 
[Data,indx]=sort(Data);

% Calculate length of Vector Input Data 
n=length(Data);

% Check if the number of Data's elements is even or odd
if mod(n,2)~=0
    % If the number is odd, split the set as explained in the oddlength
    % function
    [~,Data1,Data3,n2]=oddlength(Data,n);
    % Check if the number of the split Datas' elements is even or odd
    if mod(n2,2)~=0
        % If the number is odd, calculate the median of each split set,
        % according to the oddlength function
        Q1=oddlength(Data1,n2);
        Q3=oddlength(Data3,n2);
    else
        % If the number is even, calculate the median of each split set,
        % according to the evenlength function
        Q1=evenlength(Data1,n2);
        Q3=evenlength(Data3,n2);
    end
else
    % If the number is even, split the set as explained in the evenlength
    % function
    [~,Data1,Data3,n2]=evenlength(Data,n);
    % Check if the number of the split Datas' elements is even or odd
    if mod(n2,2)~=0
        % If the number is odd, calculate the median of each split set,
        % according to the oddlength function
        Q1=oddlength(Data1,n2);
        Q3=oddlength(Data3,n2);
    else
        % If the number is even, calculate the median of each split set,
        % according to the evenlength function
        Q1=evenlength(Data1,n2);
        Q3=evenlength(Data3,n2);
    end
end

% Calculate the value of IQR
vIQR=Q3-Q1;

% Unsort the Data
UData(indx)=Data;

% Transpose the Unsorted Data matrix
UData=transpose(UData);

% Delete the values that are outside the IQR range
idx=UData>=Q3+1.5*vIQR | UData<=Q1-1.5*vIQR;
UData(idx)=[];
Time(idx)=[];

% Save the filtered Data
FilteredData=UData;
end
function [Median,LowerSet,UpperSet,n2]=oddlength(Data,n)
%% This function calculates the Median of the Input Vector Data set, 
%% divides the set to Upper and Lower halves and calculates the length of the produced sets
%% when the number of elements is odd

% The median of an odd set is the element that divides the set into two
% sets with an equal number of elements
Median=Data((n-1)/2+1);

% The lower half of the set contains elements from the first one to the one before
% the median of the original set
LowerSet=Data(1:(n-1)/2);

% The upper half of the set contains elements from the one after
% the median to the end of the original set
UpperSet=Data((n-1)/2+2:end);

% Calculate the length of the produced sets
n2=length(LowerSet);
end
function [Median,LowerSet,UpperSet,n2]=evenlength(Data,n)
%% This function calculates the Median of the Input Vector Data set, 
%% divides the set to Upper and Lower halves and calculates the length of the produced sets
%% when the number of elements is even

% The median of an even set is the mean value of the two central elements
% that split the set to two sets with an equal number of elements
Median=(Data(n/2)+Data(n/2+1))/2;

% The lower half of the set contains elements from the first one until the
% former central element of the original set
LowerSet=Data(1:n/2);

% The upper half of the set contains elements from the first latter central element 
% until the end of the original set
UpperSet=Data(n/2+1:end);

% Calculate the length of the produced sets
n2=length(UpperSet);
end

%% 
% *To plot data filtered and unfiltered:*

function plot_differences(Data2024,fieldnames,units,data_filt,time_filt,i)
    % Plot unfiltered data
    figure;
    subplot(2,1,1); % Plot on the first row, first two columns
    plot(Data2024.Time, Data2024.(fieldnames{i}),".");
    title(fieldnames{i}); % Set title as column name
    xlabel('Time [s]'); % X-axis label
    ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label
    
    % Plot unfiltered data
    subplot(2,1,2); % Plot on the second row, first column
    plot(time_filt.(fieldnames{i}),data_filt.(fieldnames{i}),".");
    title(['Filtered ' fieldnames{i}]);
    xlabel('Time [s]'); % X-axis label
    ylabel(sprintf('%s %s', fieldnames{i}, units{i})); % Concatenate field name with unit and set it as Y-axis label
end
%% 
% 

function [] = pred_and_plot(x,y,net,i,j,fieldnames)
    preds = predict(net,x);
    figure
    x=x(:,i);
    scatter(x,y,'r')
    hold on
    hold on
    scatter(x,preds,'g')
    title('Test vs Net prediction')
    ylabel(sprintf('%s', fieldnames{j}));
    xlabel(sprintf('%s', fieldnames{i}));
    legend('Test','Net Prediction')
    grid on
end
%% 
% 

function [] = pred_and_plot2(x,y,net,i,j,fieldnames,units,center,scale)
    preds = predict(net,x);
    figure
    x = x(:,i);
    x = x*scale{:,i}+center{:,i};
    y = y*scale{:,j}+center{:,j};
    preds = preds*scale{:,j}+center{:,j};
    scatter(x,y,'r')
    hold on
    hold on
    scatter(x,preds,'g')
    title('Test vs Net prediction')
    ylabel(sprintf('%s %s', fieldnames{j},units{j}));
    xlabel(sprintf('%s %s', fieldnames{i},units{j}));
    legend('Test','Net Prediction')
    grid on
end