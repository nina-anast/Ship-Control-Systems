function [] = pred_and_plot2(normalizedTesting,net,i,j,fieldnames,units,centerValueTraining,scaleValueTraining)

%% This function predicts and plots the requested engine parameter in respect to time 
%% using testing data and trained network.
% Input the normalized testing data (normalizedTesting), the trained net (net), the requested parameter's 
% column number i, the parameter's (that was used for training) column number j, cell containing the names
% of all engine parameter's (fieldnames), cell containing the unitsof allengine parameter's (units),
% center and scale values for later de-normalization of the predicted data
% (centerValueTraining and scaleValueTraining)

% Isolate time from testing data used in plot
x = table2array(normalizedTesting(:,1));
% Isolate requested parameter from testing data
y = table2array(normalizedTesting(:,i));
% Isolate parameter used in training to predict the requested parameter
z = table2array(normalizedTesting(:,j));
% Predict the requested parameter
preds = predict(net,z);
figure;
% De-normalize time, test's requested parameter and prediction of requested
% parameter
x = x*scaleValueTraining{:,1}+centerValueTraining{:,1};
y = y*scaleValueTraining{:,i}+centerValueTraining{:,i};
preds = preds*scaleValueTraining{:,i}+centerValueTraining{:,i};
% Plot test's requested parameter and prediction of requested in respect to time
plot(x,y,'r.',x,preds,'g.')
title('Test vs Net prediction')
ylabel(sprintf('%s %s', fieldnames{i},units{i}));
xlabel(sprintf('%s %s', fieldnames{1},units{1}));
legend('Test','Net Prediction')
grid on
end

