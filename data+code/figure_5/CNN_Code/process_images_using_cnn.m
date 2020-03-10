function status = process_images_using_cnn(input_directory, ext, output_directory, ...
    cnn_net, neighborhood_size, scale, overwrite)

    status = false;

    % Fantastic, first, check to see if the output directory exists
    % Also check to see if classified_images directory exists
    
    directory_test({output_directory,...
        fullfile(output_directory, 'classified_mats')});

    % Now, get a listing of all files to process
    
    to_process = dir(fullfile(input_directory, ['*.', ext]));
    
    % Number of files
    
    num_to_process = length(to_process);
    
    % Create a cell array that will contain the filename and then the
    % tabulated data
    
    cell_collector = cell(num_to_process, 2);
    
    % Let's loop through -> this is going to take a while!
    tic
    for x = 1:num_to_process
        this_file = fullfile(input_directory, to_process(x).name);
        [~, name_only, ~] = fileparts(this_file);
        output_name = fullfile(output_directory, 'classified_mats', [name_only, '_classified.mat']);
        % If file exists and overwrite is false, then move on
        if exist(output_name, 'file') == 2 && ~overwrite
            continue
        end
        % Otherwise, let's classify!
        classified = classify_pixels(this_file,...
            scale, neighborhood_size, cnn_net);
        % Tabulate!
        freq = tabulate(classified(:));
        % Update cell_collector 
        cell_collector{x, 1} = name_only;
        cell_collector{x, 2} = freq;
        % Save the mats
        save(output_name, 'classified');
        save(fullfile(output_directory, 'stats.mat'), 'cell_collector');
        disp(['Image ', name_only, ' just processed']);
    end
    toc
    status = true;
end