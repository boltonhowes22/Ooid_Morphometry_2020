function [imgs, labels, validation_imgs, validation_labels] = produce_training_data_for_cnn(img_location,...
    input_idxs,...
    input_labels,...
    neighborhood, ...
    number_validation)
    % This function loads the image of interest and produces 4D arrays for
    % a CNN classifier.
    % neighboorhood is a 1x2 matrix that corresponds to the y x dimensions
    % of the nieghborhood around a pixel
    % Start by making sure that the neighborhood values are odd (so as to
    % center on the pixel of interest)
    %% A S S E R T I O N S
    assert(is_odd(neighborhood(1)) && is_odd(neighborhood(2)),...
        'Error! Neighborhood dimensions must be odd!');
    %% P R E - P R O C E S S I N G
    % Begin by loading in the image
    img = imread(img_location);
    % Convert to double
    img = im2double(img);
    % Get the size of the original image
    [img_y, img_x, channels] = size(img);
    % Next, pad the image by the neighborhood size
    padded_img = padarray(img, neighborhood, 'symmetric');
    %% SPLIT DATA
    % Here, we split the input data into input and validation data
    val_idx = randperm(size(input_idxs,1), number_validation);
    %% N E I G H B O R H O O D S
    %% X and Y half size for neighborhoods
    y_half = (neighborhood(1) - 1) / 2;
    x_half = (neighborhood(2) - 1) / 2;
    % Number of input idxs
    number_idxs = length(input_idxs);
    imgs = zeros(neighborhood(1), neighborhood(2), channels, number_idxs - number_validation);
    labels = cell(number_idxs - number_validation, 1);
    validation_imgs = zeros(neighborhood(1), neighborhood(2), channels, number_validation);
    validation_labels = cell(number_validation, 1);
    img_counter = 1;
    val_counter = 1;
    % Now, for each input_idx
    for x = 1:number_idxs
        % Convert the original index into a subscript so that we can convert
        % it into a padded array index. Note that we only need to do this
        % on a single 2D channel.
        [idx_row, idx_col] = ind2sub([img_y, img_x], input_idxs(x));
        % Convert to a padded_img x and y
        padded_idx_row = idx_row + neighborhood(1);
        padded_idx_col = idx_col + neighborhood(2);
        % Declare data matrix 
        data = zeros(neighborhood(1), neighborhood(2), channels);
        % For each channel, add to data
        for y = 1:channels
            data(:,:,y) = padded_img(padded_idx_row - y_half : padded_idx_row + y_half, ...
                padded_idx_col - x_half : padded_idx_col + x_half, y);
        end
        % Checking to see if this is a validation image
        if ismember(x, val_idx)
            validation_imgs(:,:,:,val_counter) = data;
            validation_labels{val_counter} = input_labels{x};
            val_counter = val_counter + 1;
        else
            imgs(:,:,:,img_counter) = data;
            labels{img_counter} = input_labels{x};
            img_counter = img_counter + 1;
        end
        data = [];
    end
end

function answer = is_odd(val)
    % Function to test if a value is odd. Ridiculous that I'm writing it, I
    % know.
    if (mod(val, 2))
        answer = true;
    else
        answer = false;
    end
end