function classified_image = classify_pixels(image, scale, neighborhood, neural_net)
    % This image accepts a class list, an image, and a CNN network and
    % processes each pixel of the image, assigning one class to each pixel.
    % The output is a label image where the numbers correspond to the
    % class_list index of the assigned class. This function assumes that
    % the neighborhood and class_list inputs match those provided in the
    % training set.
    
    % First, load the image - REMOVE THIS UINT 8 CONVERSION!
    img = imread(image);
    
    % Rescale if scale is not 1
    
    if scale ~= 1
        img = imresize(img, scale);
    end
     
    % Convert to double
    % img = im2double(img);
    
    % Convert to uint8?
    % img = uint8(img);

    % Get the size of the original image
    
    [img_row, img_col, ~] = size(img);
    
    % Neighborhood vals
    
    neighborhood_row = neighborhood(1);
    neighborhood_col = neighborhood(2);
   
    % Calculate the neighborhood half sizes
    row_half = (neighborhood_row - 1) / 2;
    
    col_half = (neighborhood_col - 1) / 2;

    % Then, pad the image
    
    padded_img = padarray(img, [row_half, col_half], 'symmetric');
    
    [padded_row, padded_col, img_channel] = size(padded_img);
    
    % Square subimage size
    
    sub_image_dim = 500;
    
    % Next, determine how many subimages you will need in the y and the x
    
    row_sub_image = ceil(img_row / sub_image_dim);
    
    col_sub_image = ceil(img_col / sub_image_dim);
    
    % Now, let's produce a start and end matrix for row and col
    
    % Begin by producing start and end matricies for the unpadded image
    
    start_row = (((1:row_sub_image) - 1) * sub_image_dim) + 1;

    end_row = min(start_row + sub_image_dim - 1, img_row);
    
    start_col = (((1:col_sub_image) - 1) * sub_image_dim) + 1;
    
    end_col = min(start_col + sub_image_dim - 1, img_col);

    % Now, shift the end values by 2 * half neighborhood distance -> this
    % is because the start values would have to shift by start_value +
    % half_neighborhood - half_neighborhood, so effectively by 0. End
    % values shift by end_value + half_neighborhood + half_neighborhood

    end_row = end_row + (2 * row_half);

    end_col = end_col + (2 * col_half);
  
    % Empty classified image
    
    classified_image = zeros(img_row, img_col);
    
    tic

    for this_row = 1:length(start_row)
        for this_col = 1:length(start_col)
            classified_vals = semanticseg(...
                padded_img(start_row(this_row):end_row(this_row), ...
                start_col(this_col):end_col(this_col),...
                :), ...
                neural_net, ...
                'OutputType', 'double');
            classified_image(((this_row - 1) * sub_image_dim) + 1 : min(this_row * sub_image_dim, img_row), ...
                ((this_col - 1) * sub_image_dim) + 1 : min(this_col * sub_image_dim, img_col)) = classified_vals;
        end
    end
    
    toc

end