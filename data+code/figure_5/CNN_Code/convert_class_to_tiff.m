function convert_class_to_tiff(input_dir, output_dir)
    % Check to make sure that output_dir exists
    directory_test({output_dir});
    % Get list of all mat files
    to_process = dir(fullfile(input_dir, '*.mat'));
    % Lets goooo
    parfor x = 1:length(to_process)
        % Filepath
        this_path = fullfile(to_process(x).folder, to_process(x).name);
        % This name only
        [~, name_only, ~] = fileparts(this_path);
        % Load the file
        loaded_struct = load(this_path);
        % Get fieldnames
        fnames = fieldnames(loaded_struct);
        % Set data
        data = loaded_struct.(fnames{1});
        % Now, write
        imwrite(uint8(data), fullfile(output_dir, [name_only, '.tif']));
    end
end