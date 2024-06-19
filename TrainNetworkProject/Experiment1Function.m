function [layers,options] = Experiment1Function(params)
    switch params.layer
        case "basic"
            layers = [
                featureInputLayer(3)
                fullyConnectedLayer(20)
                reluLayer
                dropoutLayer(0.05)
                fullyConnectedLayer(1)
                regressionLayer];
        case "2-layer"
            layers = [
                featureInputLayer(3)
                fullyConnectedLayer(20)
                tanhLayer
                fullyConnectedLayer(1)
                regressionLayer];
        otherwise
            error("Unknown network type.");
    end
    options = trainingOptions('adam',...
    'Shuffle','every-epoch',...
    'MiniBatchSize',1000, ...
    'MaxEpochs',10, ...
    'InitialLearnRate',0.0001, ...
    'Plots','training-progress', ...
    'Verbose',false);
end