function produce_training_data_for_cnn_v2(img_location,...
    input_idxs,...
    output_labels,...
    neighborhood, ...
    convert_to_double, ...
    output_directory, ...
    output_ext)
    % This function loads the image of interest and produces 4D arrays for
    % a CNN classifier.
    % neighboorhood is a 1x2 matrix that corresponds to the y x dimensions
    % of the nieghborhood around a pixel
    % Start by making sure that the neighborhood values are odd (so as to
    % center on the pixel of interest)
    %% A S S E R T I O N S
    % Check if neighborhood size is odd. If not, throw an error. 
    assert(is_odd(neighborhood(1)) && is_odd(neighborhood(2)),...
        'Error! Neighborhood dimensions must be odd!');
    % What are the unique labels?
    unique_labels = unique(output_labels);
    % The count
    num_unique_labels = length(unique_labels);
    % What are the label directories?
    label_directories = cellfun(@(x) fullfile(output_directory, x), unique_labels, 'UniformOutput', false);
    % Do said directories (including the overarching output_directory) exist?
    directory_test(vertcat(output_directory, label_directories));
    % Next, we'd like to set up image names correctly
    labels_number = struct;
    for label_ct = 1:num_unique_labels
        % Add this label to labels_number
        labels_number.(unique_labels{label_ct}) = struct;
        labels_number.(unique_labels{label_ct}) = 1;
        % Okay, let's get a directory listing 
        list_files = dir(fullfile(output_directory, unique_labels{label_ct}, ...
            ['*', output_ext]));
        if ~isempty(list_files)
            % Sort these
            file_names = {list_files.name};
            [~, sorted_idx] = sort(file_names);
            sorted = list_files(sorted_idx);
            % Get last value
            last_value = sorted(end).name;
            % Break into parts
            [~, last_name, ~] = fileparts(last_value);
            % Now, set this label_number to the last_name + 1;
            labels_number.(unique_labels{label_ct}) = str2num(last_name) + 1;
        end
    end
    %% P R E - P R O C E S S I N G
    % Begin by loading in the image
    img = imread(img_location);
    % Convert to double?
    % Note that, for jpg inputs, this conversion is important...
    % Look into why...
    if convert_to_double 
        img = im2double(img);
    end
    % Get the size of the original image
    [img_y, img_x, channels] = size(img);
    % Next, pad the image by the neighborhood size
    padded_img = padarray(img, neighborhood(1:2), 'symmetric');
    %% N E I G H B O R H O O D S
    %% X and Y half size for neighborhoods
    y_half = (neighborhood(1) - 1) / 2;
    x_half = (neighborhood(2) - 1) / 2;
    % Number of input idxs
    number_idxs = length(input_idxs);
    progress_bar = waitbar(0, 'Processing Images...');
    tic 
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
        % Note that we'll need to figure out how to correctly deal with
        % different types of data. Right now, don't declare as uint16
        if convert_to_double
            data = zeros(neighborhood(1), neighborhood(2), channels);
        else
            data = zeros(neighborhood(1), neighborhood(2), channels, 'uint16');
        end
        % For each channel, add to data
        for y = 1:channels
            data(:,:,y) = padded_img(padded_idx_row - y_half : padded_idx_row + y_half, ...
                padded_idx_col - x_half : padded_idx_col + x_half, y);
        end
        % This label
        this_label = output_labels{x};
        % What is this label number?
        this_number = labels_number.(this_label);
        % Write this to the output_directory corresponding to this label
        if convert_to_double
            imwrite(data, fullfile(output_directory, this_label, [sprintf('%06d', this_number), output_ext]));
        else
            write_tiff(data, fullfile(output_directory, this_label, [sprintf('%06d', this_number), output_ext]));
        end
        % Update the labels_number for this label
        labels_number.(this_label) = this_number + 1;
        % Update a progress bar
        waitbar(x/number_idxs, progress_bar, ['Image ', num2str(x), ' of ', num2str(number_idxs), ' written']);
    end
    close(progress_bar);
    toc
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