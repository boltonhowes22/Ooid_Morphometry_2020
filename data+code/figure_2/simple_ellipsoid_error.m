function [major_length, minor_length] = simple_ellipsoid_error(aspect_ratio_ac, aspect_ratio_bc , max_num_layers, number_obs, nuc_bool, nuc_percent)
% Makes a defined elliposoid, then cuts with a randomly oriented and placed
% plane to observe the errors of 2D measurements on 3D objects. 
% a<b<c
% IN
% aspect_ratio_ac: ratio of a to c
% aspect_ratio_bc: ratio of b to c
% number_obs: How many measurements per rep
% max_num_layers: How many growth bands are in the ooid?
% nuc_bool: boolean to indicate presence or absence of a nucleus
% nuc_percent: percent of entire ellipsoid occupied by a nucleus
%
% OUT
% major_length: the measured major axis lengths
% minor_length: the measured minor axis lengths
%
% Bolton Howes
% January 2020
%% Dimension

% Make outer ooid
long_axis = 50; 
intermediate_axis = long_axis * aspect_ratio_bc;
short_axis = long_axis * aspect_ratio_ac;

% Make inner ooid
if nuc_bool == 1
    long_axis2 = 50 * nuc_percent; 
    intermediate_axis2 = long_axis2 * aspect_ratio_bc;
    short_axis2 = long_axis2 * aspect_ratio_ac;
end
%%
% preallocate size of storage
if max_num_layers ~= 1
%     major_length = zeros(length(1), reps, 2 - 1);
%     minor_length = zeros(length(1), reps, 2 - 1);
%     eccentricity = zeros(length(1), reps, 2 - 1);
else
    major_length = zeros(number_obs,1);
    minor_length = zeros(number_obs,1);
    %eccentricity = zeros(number_obs, reps);
end
%% Create meshgrid
A = -long_axis-15:long_axis+15;
[X,Y,Z] = meshgrid(A,A,A);

%% In the ellipsoid

if max_num_layers ~=1
    in_ellipsoid1 = (X./long_axis).^2 + (Y./intermediate_axis).^2 + (Z./short_axis).^2 <= 1;
    in_ellipsoid2 = (X./long_axis2).^2 + (Y./intermediate_axis2).^2 + (Z./short_axis2).^2 <= 1;
    in_ellipsoid2 = in_ellipsoid2 .* 100;
    in_ellipsoid1 = in_ellipsoid1 .* 100;
    in_ellipsoid = (in_ellipsoid1 + in_ellipsoid2)./255;
else 
    in_ellipsoid = (X./long_axis).^2 + (Y./intermediate_axis).^2 + (Z./short_axis).^2 <= 1;
end



%% Slice the ooids
sliced_imgs = shape_slicer(in_ellipsoid, number_obs);

%% Count number of layers in the ooid



tracker_complete_layers = 0;

if max_num_layers ~=1

    num_layers = zeros(number_obs, 1);
    if number_obs > 10

    end

    for n = 1:length(sliced_imgs)
        current_im = sliced_imgs{n};

        % Histogram of pixel values
        p = imhist(current_im);

        % Find peaks on the histogram
        [pks, ~] = findpeaks(p, 'MinPeakHeight', 1000);

        % Number of classes in the kmeans is the number of peaks in the histogram + 1 for the
        % background 
        num_classes = length(pks) + 1;

        % How many layers are being collected? This will be used to test the
        % likelihood of intersecting the nucleus of an ooid
        num_layers(n) = length(pks);

    %     % Convert image to uint8 for the k means clustering
    %     d = linspace(min(current_im(:)), max(current_im(:)),256);
    %     im1 = uint8(arrayfun(@(x) find(abs(d(:)-x) == min(abs(d(:)-x))),current_im));
    %     
    %     % K means
    %     pixel_labels = imsegkmeans(im1,num_classes);
    %     
        % Visualize the k means
        %B = labeloverlay(current_im, pixel_labels); 
        %figure(9)
        %imshow(B)

        % Pull all samples that intersect all layers
        if num_layers(n) ==  max_num_layers
            tracker_complete_layers = 1 + tracker_complete_layers;
            all_layers{tracker_complete_layers} = sliced_imgs{n};
        end


        if number_obs > 10

        end

    end
    if number_obs > 10
        %close(h) 
    end

else
    % There is only one layer, and this step is unecessary
end

%% Do the measuring 

% Keep track
count1 = 0;

% Begin measuring


if max_num_layers ~= 1
    for j = 1:length(all_layers)

        % Get Current Image
        current_im = all_layers{j};

        % Convert image to uint8
        d = linspace(min(current_im(:)), max(current_im(:)),256);
        im1 = uint8(arrayfun(@(x) find(abs(d(:)-x) == min(abs(d(:)-x))),current_im));

        % K means
        pixel_labels = imsegkmeans(im1,num_classes);

        for n = 2:num_classes


            % Mask the pixel lables
            mask = pixel_labels == n;
            %im2 = im1 < 5 == 0;
            growth_band =  uint8(im1) .* uint8(mask);
            growth_band = imbinarize(growth_band);

            % Make Get the convex image to measure eccentricity 
            convex_image = regionprops(growth_band, 'ConvexImage');
            
            % Pad the image with 0's to make the measurement easier
            %convex_image.ConvexImage = padarray(convex_image.ConvexImage, [5,5]);
            
            
            if length(convex_image) ~=1
                keyboard
            end
            % Check to make sure the correct number of peaks were found 
            if length(convex_image) == 1

                
                
                 % Get boundaries
                [boundary, ~] = bwboundaries(convex_image.ConvexImage);

                coordinates = [boundary{1}(:,2), boundary{1}(:,1)];
                coord_error_count = 0;
                if sum(coordinates(:)) == 4
                    disp('skip coordinate error')
                    coord_error_count = coord_error_count + 1;
                    continue
                end
                % Rotate to orient along major axis
                rotated = pca_rotate(coordinates);

                minor_length_temp = abs(min(rotated(:,1))) +  max(rotated(:,1));

                major_length_temp = abs(min(rotated(:,2))) +  max(rotated(:,2));

   
                
                if major_length_temp > 110
                    count1 = count1 +1; 
                    disp(['oops 1 #' num2str(count1)])

                else

                    % Record Values
                    major_length(j, n-1) = major_length_temp;
                    minor_length(j, n-1) = minor_length_temp;

                end


            else

                % Sometimes the wrong number of peaks are found 
                count1 = count1 +1; 
                disp(['oops #' num2str(count1)])
            end


        end

    end

else
    % There is only one layer, so the K-Means isn't necessary

    for j=1:number_obs

        % Make Get the convex image to measure eccentricity 
        convex_image = regionprops(sliced_imgs{j}, 'ConvexImage', 'Eccentricity');

        % Produce Convex Image and Measure Eccentricity
        % convex_image_pixel_list = regionprops('table',convex_image.ConvexImage, 'PixelList');

        % Get boundaries
        [boundary, ~] = bwboundaries(convex_image.ConvexImage);
        
        coordinates = [boundary{1}(:,2), boundary{1}(:,1)];
        coord_error_count = 0;
        if sum(coordinates(:)) == 4
            disp('skip coordinate error')
            coord_error_count = coord_error_count + 1;
            continue
        end
        % Rotate to orient along major axis
        rotated = pca_rotate(coordinates);

        minor_length(j) = abs(min(rotated(:,1))) +  max(rotated(:,1));

        major_length(j) = abs(min(rotated(:,2))) +  max(rotated(:,2));



        %{
        try
            % Measure Feret Dimensions
            feret_table = feretProperties(convex_image_pixel_list);
                        % Record Values
            major_length(j, rep) = feret_table.MaxFeretDiameter;
            minor_length(j, rep) = feret_table.MinFeretDiameter;
            %eccentricity(j, rep) = convex_image_pixel_list.Eccentricity;

        catch ME
            keyboard
            disp('Pesky Error was skipped!')
            continue
        end
        %}
    end


end


end
