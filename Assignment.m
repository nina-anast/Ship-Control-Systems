%% *Εργασία Ειδικά Συστήματα Πλοίου*
% 
%% *Ερώτημα 1:* Ανάλυση των δεδομένων με τη χρήση διαφόρων γραφικών παραστάσεων
% 

clear all
warning("off")
% load all data
Data2024 = readtable("Data2024.csv",'TrailingDelimitersRule','ignore')
% Extract the units
% Open the file for reading
fid = fopen('Data2024.csv', 'r');
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
clearvars num_cols num_rows num_plots i ans
% Ανάλυση Δεδομένων

i = 2;
figure; % Create a new figure
plot(data.Time, data.(fieldnames{i}),".");
title(fieldnames{i}); % Set title as column name
xlabel('Time [s]'); % X-axis label

% Concatenate field name with unit and set it as Y-axis label
ylabel(sprintf('%s %s', fieldnames{i}, units{i}));
clearvars i
%% 
% *NOx:* Έχουν περιοδική κίνηση με απόσβεση, φθίνουσα. Περίοδος είναι 225s. 
% Η μέση τιμή μικραίνει με την πάροδο του χρόνου, ξεκινάει από 600 ppm και καταλήγει 
% εως και 200 ppm.
% 
% *Fuel Consumption:* Έχουν περιοδική κίνηση με αυξανόμενο πλάτος (ίδια περίοδος 
% με ΝΟx 180-405s η πρώτη περίοδος). Οι μέγιστες τιμές ξεκινάνε από 20kg/h και 
% φτάνουν ως 65kg/h.
% 
% *Λάμδα:* Έχει τιμές μικρότερες του 20 εως την στιγμή 2619s που γίνεται το 
% πρώτο peak. Έχει συνολικά 8 peaks με τα μέγιστα να ξεπερνάνε την τιμή 130.
%% Ερώτημα 2: Φιλτράρισμα και μετασχηματισμός των δεδομένων.
% Το παρακάτω μπορεί να χρησιμοποιηθεί για να παράγουμε τα φιλτραρισμένα δεδομένα, 
% επιλέγοντας κάθε φορά το σωστό threshold.

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
% Αποθήκευση Φιλτραρισμένων Δεδομένων
% Για να αναπαράγουμε τα παραπάνω δεδομένα θα αποθηκεύσουμε τα thresholds που 
% παράγαμε και θα μπορούμε να φορτώνουμε ξανά τα παραπάνω δεδομένα.

% writematrix(thresholds,"thresholds.csv");
% clearvars ans 
% Ανάκτηση Δεδομένων
% Για να ξαναποκτήσουμε τα φιλτραρισμένα δεδομένα, μπορούμε να τρέξουμε τον 
% παρακάτω κώδικα:

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
% Μέθοδος *Interquartile Range*:
% 
% Apply to:
%% 
% * Fuel Consumption (3)
% * Intake Pressure (6)
% * Torque Reference (7)
% * Rot Speed (8)
% * Engine Torque (9)
% * (EGR Command(10) ->ιδιο με threshold)
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
% *Τελικό φιλτράρισμα:*

% Initialize a struct to store filtered values
data_filt = struct();
time_filt = struct();

% Read thresholds from the CSV file
thresholds = readmatrix('thresholds.csv');

for j = 2:numel(fieldnames)
    if j == 2 || j == 4 || j == 10
        [filteredData, filteredTime, ~] = chauvenetFilter2(Data2024, fieldnames, units, j, thresholds(j));
        data_filt.(fieldnames{j}) = filteredData;
        time_filt.(fieldnames{j}) = filteredTime;
    elseif j == 3 || j == 6 || j == 7
        [filteredData,filteredTime] = IQR(Data2024.(fieldnames{j}),Data2024.Time);
        data_filt.(fieldnames{j}) = filteredData;
        time_filt.(fieldnames{j}) = filteredTime;
    elseif j == 5 || j == 8 || j == 9 || j == 11
        % Filter the data using rmoutliers
        [data_filt.(fieldnames{j}), idx] = rmoutliers(Data2024.(fieldnames{j}), "median");
        
        % Retain corresponding times using logical indexing
        time_filt.(fieldnames{j}) = Data2024.Time;
        time_filt.(fieldnames{j})(idx) = [];
    end
end

clearvars filteredData filteredTime j idx
%%
plot_differences(Data2024,fieldnames,units,data_filt,time_filt,2)
%% 
% Για την επιλογή των *thresholds*:

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
% Για την *ανάκτηση* των φιλτραρισμένων δεδομένων:

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
% Μέθοδος *Interquartile Range*:
% 
% <https://www.indeed.com/career-advice/career-development/outliers-statistics 
% https://www.indeed.com/career-advice/career-development/outliers-statistics>

function [FilteredData,Time]=IQR(Data,Time)
%% Remove outliers from vectors with interquartile range method

% Sort Input Data in ascending order 
[Data,indx]=sort(Data);

% Calculate length of Vector Input Data 
n=length(Data);

% Check if the number of Data points is odd or even
if mod(n,2)~=0
    % Number of Data is odd
    % Divide the Data into two sets. The lower half of the data set saves
    % the values until the element before the median, while the upper
    % half saves values starting from the element after the median until
    % the end. The median is the element that divides the Data set into
    % two separate sets with the same number of elements.

    % Define lowet half of Data set
    Data1=Data(1:(n-1)/2);
    % Calculate length of lower half set
    n1=length(Data1);
    % Define upper half set of Data
    Data3=Data((n-1)/2+2:end);
    % Calculate length of upper half set
    n3=length(Data3);
    
    % The corresponding sets have an even number of elements
    % The median of an even set is the mean of the two central elements,
    % which divide the set in two separate sets with the same number of
    % elements

    % Calculate the first quartile (the median of the lower half of the data set)
    Q1=(Data1(n1/2)+Data1(n1/2+1))/2;
    % Calculate the third quartile (the median of the upper half of the data set)
    Q3=(Data3(n3/2)+Data3(n3/2+1))/2;
else
    % Number of Data is even
    % Divide the Data into two sets with the same number of elements

    % Define lowet half of Data set
    Data1=Data(1:n/2);
    % Calculate length of lower half set
    n1=length(Data1);
    % Define upper half set of Data
    Data3=Data(n/2+1:end);
    % Calculate length of upper half set
    n3=length(Data3);

    % The corresponding sets have an even number of elements

    % Calculate the first quartile (the median of the lower half of the data set)
    Q1=(Data1(n1/2)+Data1(n1/2+1))/2;
    % Calculate the third quartile (the median of the upper half of the data set)
    Q3=(Data3(n3/2)+Data3(n3/2+1))/2;
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