function [cnn_net, training_info] = train_cnn_multiple_images(...
    training_directory, neighborhood_size, ...
    training_image_list, training_ext, convert_to_double, imds_storage_directory, ...
    generate_imds, label_map, network_directory)
    % If generate_imds is set to true, we're going to create the imds
    % directories and associated training data
    if generate_imds
        for this_file = 1:length(training_image_list)
            % This image 
            this_image = fullfile(training_directory, training_image_list{this_file});
            % Take the image and add _training to the end of the name
            [this_dir, this_name, ~] = fileparts(this_image);
            % Where is the training image?
            training_image_path = fullfile(this_dir, [this_name, '_training', training_ext]);
            training_image = imread(training_image_path);
            % Turn your labelled image into the input_idxs and the
            % output_labels for futher processing
            % The third argument is the number of training pixels per class
            [input_idxs, output_labels] = hand_class_to_idx(training_image, label_map, 8.5e4);
            % Here, we want to create a imds for further training. Now, one
            % thing to think about is how long it takes to write (and read!)
            % that many training data to file.
            produce_training_data_for_cnn_v2(this_image, ...
                input_idxs, output_labels, neighborhood_size, ...
                convert_to_double, imds_storage_directory, '.tif');
        end
    end
    
    % Else, let's go ahead and declare the imds_storage_directoy as an imds
    data_store = imageDatastore(imds_storage_directory, 'LabelSource', ...
        'foldernames', 'IncludeSubfolders', true);
    % Next, let's split the data into training and validation
    [training, validation] = splitEachLabel(data_store, 0.75, 0.25, 'randomize');
    
    % Awesome, let's create, and then train, the neural network...
    
    % Begin with an input layer
    input_layer = imageInputLayer(neighborhood_size, 'Name', 'in');
    
    % Local Convolution Layers
    local_convolution_layers = [convolution2dLayer([3, 3], 64, 'Name', 'conv_1');
        reluLayer('Name', 'relu_1');
        maxPooling2dLayer([3, 3], 'Name', 'maxpool_1');
        convolution2dLayer([2, 2], 64, 'Name', 'conv_2');
        reluLayer('Name', 'relu_2');
        maxPooling2dLayer([2, 2], 'Name', 'maxpool_2');];
    
    % Global Convolution Layers
    global_convolution_layers = [
        convolution2dLayer([7, 7], 160, 'Name', 'conv_3'); 
        reluLayer('Name', 'relu_3');];

    % Final Convolution Layers
    final_convolution_layers = [depthConcatenationLayer(2, 'Name', 'depth_1'); 
        convolution2dLayer([5, 5], length(label_map{1}), 'Name', 'conv_4');
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
    
    % Now, let's set up the training parameters
    
    % Check to make sure that the network directory exists!
    
    directory_test({network_directory});
   
    miniBatchSize = 32;
    
    numValidationsPerEpoch = 10;
    
    validationFrequency = floor(size(training.Labels, 1)/miniBatchSize/numValidationsPerEpoch);
    
    options = trainingOptions('sgdm',...
        'LearnRateSchedule', 'piecewise',...
        'InitialLearnRate', 1e-6,...
        'MaxEpochs', 25,...
        'ValidationPatience', 25,...
        'MiniBatchSize', miniBatchSize,...
        'ValidationData', validation,...
        'ValidationFrequency', validationFrequency,...
        'Plots', 'training-progress',...
        'CheckpointPath', network_directory);
    
    
    % Now, train
    
    [cnn_net, training_info] = trainNetwork(training, layer_graph, options);

end