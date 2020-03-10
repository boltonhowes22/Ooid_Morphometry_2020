function nn_controller(image_directory, output_directory)
    % Create directories if they don't exist
    if ~exist(output_directory, 'dir')
        mkdir(output_directory);
    end
    label_matrix_dir = fullfile(output_directory, 'labelled_matrix');
    if ~exist(label_matrix_dir, 'dir')
        mkdir(label_matrix_dir);
    end
    class_dir = fullfile(output_directory, 'classified_matrix');
    if ~exist(class_dir, 'dir')
        mkdir(class_dir);
    end
    output_image_dir = fullfile(output_directory, 'output_matrix');
    if ~exist(output_image_dir, 'dir')
        mkdir(output_image_dir);
    end
    % First, read through the image directory and list all the tiffs
    img_list = dir(fullfile(image_directory, '*.tif'));
    parfor x = 1:length(img_list)
        fprintf(['Running iteration', num2str(x), '\n']);
        % Define the image name
        [~, name_only, ~] = fileparts(img_list(x).name);
        [label_matrix, input_matrix] = create_superpixel(fullfile(image_directory, img_list(x).name));
        classified_matrix = nntest(input_matrix);
        output_image = replace_regions(classified_matrix, label_matrix);
        parsave(label_matrix_dir, [name_only, '.mat'], label_matrix);
        parsave(class_dir, [name_only, '.mat'], classified_matrix);
        % save(fullfile(label_matrix_dir, [name_only, '.mat']), 'label_matrix');
        % save(fullfile(class_dir, [name_only, '.mat']), 'classified_matrix');
        imwrite(output_image, fullfile(output_image_dir, [name_only, '_classified.tif']));
    end
end