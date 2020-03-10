function [cnn_net, training_info] = create_cnn(input_layer_size,...
    number_classes,...
    training_imgs,...
    training_labels, validation_imgs, validation_labels, save_location)
    % This function produces a CNN (convolutional neural network) for
    % per-pixel training of GIRI images
    % Note that input_layer_size is dicatated by the user (size of the
    % mxn neighborhood
    % Number of classes is self explanatory
    % training_data MUST be a imageDatastore object OR images as a 4D array
    % (where the first three dimensions correspond to height, width,
    % channel, respectively, and the fourth dimension corresponds to each
    % image).
    % Same goes for testing_data
    %% DIRECTORY TESTING
    directory_test({save_location});
    %% A R C H I T E C T U R E 
    % Here, we define the neural network architecture
    %{
    layers = [imageInputLayer(input_layer_size, 'DataAugmentation', 'randcrop');
    convolution2dLayer([11,11], 16);
    reluLayer();
    convolution2dLayer([11,11], 4);
    reluLayer();    
    % convolution2dLayer([1,1], 16);
    % reluLayer();
    fullyConnectedLayer(4);
    softmaxLayer();
    pixelClassificationLayer()];
    %}
    %{
    % Here's a per pixel one
    layers = [imageInputLayer(input_layer_size, 'DataAugmentation', 'randcrop');
        convolution2dLayer(3, 16, 'Stride', 2, 'Padding', 1);
        reluLayer();
        transposedConv2dLayer(3, 4, 'Stride', 2, 'Cropping', 1);
        reluLayer();
        fullyConnectedLayer(4);
        softmaxLayer();
        pixelClassificationLayer()];
    %}
    % Here's another one
    
    % Begin with an input layer
    
    input_layer = imageInputLayer(input_layer_size, 'Name', 'in');
    
    % Local Convolution Layers
    
    local_convolution_layers = [convolution2dLayer([7,7], 64, 'Name', 'conv_1');
        reluLayer('Name', 'relu_1');
        maxPooling2dLayer([4,4], 'Name', 'maxpool_1');
        convolution2dLayer([3,3], 64, 'Name', 'conv_2');
        reluLayer('Name', 'relu_2');
        maxPooling2dLayer([2,2], 'Name', 'maxpool_2');];
    
    % Global Convolution Layers
    
    global_convolution_layers = [
        convolution2dLayer([13,13], 160, 'Name', 'conv_3'); 
        reluLayer('Name', 'relu_3');];

    % Final Convolution Layers
    
    final_convolution_layers = [depthConcatenationLayer(2, 'Name', 'depth_1'); 
        convolution2dLayer([21,21], number_classes, 'Name', 'conv_4');
        softmaxLayer('Name', 'softmax_1');
        pixelClassificationLayer('Name', 'px');];
    
    layer_graph = layerGraph;
    
    layer_graph = addLayers(layer_graph, input_layer);
    
    layer_graph = addLayers(layer_graph, local_convolution_layers);
    
    layer_graph = addLayers(layer_graph, global_convolution_layers);
    
    layer_graph = addLayers(layer_graph, final_convolution_layers);
    
    % Connect!
    
    layer_graph = connectLayers(layer_graph, 'in', 'conv_1');
    
    layer_graph = connectLayers(layer_graph, 'in', 'conv_3');
    
    layer_graph = connectLayers(layer_graph, 'maxpool_2', 'depth_1/in1');
    
    layer_graph = connectLayers(layer_graph, 'relu_3', 'depth_1/in2');

    %{
    layers = [
        imageInputLayer(input_layer_size, 'DataAugmentation', 'randcrop');
        convolution2dLayer([7,7], 64);
        maxPooling2dLayer([4,4]);
        convolution2dLayer([3,3], 64);
        maxPooling2dLayer([2,2]);
        convolution2dLayer([21,21], 4);
        softmaxLayer();
        pixelClassificationLayer();
        ];
    %}

    % Options
    % options = trainingOptions('sgdm', 'MaxEpochs', 1,...
    %    'InitialLearnRate', .01);
    miniBatchSize = 128;
    numValidationsPerEpoch = 10;
    validationFrequency = floor(size(training_imgs,4)/miniBatchSize/numValidationsPerEpoch);
    options = trainingOptions('sgdm',...
        'LearnRateSchedule','piecewise',...
        'LearnRateDropFactor',0.2,...
        'LearnRateDropPeriod',5,...
        'InitialLearnRate', .01,...
        'MaxEpochs',25,...
        'ValidationPatience', 25,...
        'MiniBatchSize', miniBatchSize,...
        'ValidationData',{validation_imgs, validation_labels},...
        'ValidationFrequency',validationFrequency,...
        'Plots','training-progress',...
        'CheckpointPath',save_location);
    % Train the network
    [cnn_net, training_info] = trainNetwork(training_imgs, training_labels, layer_graph, options);
end