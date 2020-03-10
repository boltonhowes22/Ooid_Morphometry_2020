function [input_idxs, output_labels] = hand_class_to_idx(label_image, label_map, max_per_class)
    % We want to list input_idxs and their associated input_labels
    % Okay, for each value, find the indicies
    input_idxs = [];
    output_labels = {};
    for x = 1:length(label_map{1})
        to_find = label_map{1}(x);
        matched_idx = find(label_image == to_find);
        if max_per_class ~= Inf
            % To draw from
            num_original = size(matched_idx,1);
            % Randomly draw max_per_class values -> this ensures we don't
            % sample from just one part of the image in the case that the
            % max_per_class is less than the number of traced/marked pixels
            % of that class
            to_draw = randperm(num_original, min(num_original, max_per_class));
            matched_idx = matched_idx(to_draw);
        end
        input_idxs = [input_idxs; matched_idx];
        output_labels = [output_labels; repmat({label_map{2}{x}}, length(matched_idx), 1)];
    end
end

