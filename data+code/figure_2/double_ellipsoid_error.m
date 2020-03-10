function [major_length, minor_length] = double_ellipsoid_error(aspect_ratio, max_num_layers, number_obs, reps, nuc_bool, nuc_major)
% Variables
% aspect_ratio = .75;             % Aspect ratio defined as percent smaller than long axis
% reps = 50;                     % How many times we check the median
% number_obs = 10;               % How many measurements per rep
% max_num_layers = 1;             % How many growth bands are in the ooid?
%% Dimension

% Make outer ooid
long_axis = 50; 
medium_axis = long_axis * aspect_ratio; 

% Make inner ooid
long_axis2 = nuc_major; 
medium_axis2 = long_axis2 * aspect_ratio; 


%%
% preallocate size of storage

major_length = zeros(number_obs, reps);
minor_length = zeros(number_obs, reps);


%% Create meshgrid
A = -long_axis-15:long_axis+15;
[X,Y,Z] = meshgrid(A,A,A);

%% In the ellipsoid


in_ellipsoid1 = (X./long_axis).^2 + (Y./medium_axis).^2 + (Z./medium_axis).^2 <= 1;
in_ellipsoid2 = (X./long_axis2).^2 + (Y./medium_axis2).^2 + (Z./medium_axis2).^2 <= 1;
in_ellipsoid2 = in_ellipsoid2 .* 100;
in_ellipsoid1 = in_ellipsoid1 .* 100;
in_ellipsoid = (in_ellipsoid1 + in_ellipsoid2)./255;




good_measurements = 0;

for rep = 1:reps
    tracker_complete_layers = 0;
    
    while tracker_complete_layers < number_obs
    %% Slice the ooids
    sliced_imgs = shape_slicer(in_ellipsoid, number_obs);

    %% Count number of layers in the ooid


        

        if max_num_layers == 2

            num_layers = zeros(number_obs, 1);


            for n = 1:length(sliced_imgs)
                current_im = sliced_imgs{n};

                % Histogram of pixel values
                p = imhist(current_im);

                % Find peaks on the histogram
                [pks, ~] = findpeaks(p, 'MinPeakHeight', 1);

                % Number of classes in the kmeans is the number of peaks in the histogram + 1 for the
                % background 
                num_classes = length(pks) + 1;

                % How many layers are being collected? This will be used to test the
                % likelihood of intersecting the nucleus of an ooid
                num_layers(n) = length(pks);

                % Pull all samples that intersect all layers
                if num_layers(n) ==  max_num_layers
                    tracker_complete_layers = 1 + tracker_complete_layers;
                    all_layers{tracker_complete_layers} = sliced_imgs{n};
                else
                    continue
                end
            end


        else
            % There is only one layer, and this step is unecessary
        end
    end
all_layers =  all_layers(1:number_obs);

    %% Do the measuring 

    % Keep track
    count1 = 0;

    % Begin measuring

    
    if max_num_layers ~= 1
        for j = 1:length(all_layers)

            % Get Current Image-- All of these will have both layers!
            current_im = all_layers{j};

            % Convert image to uint8
            d = linspace(min(current_im(:)), max(current_im(:)),256);
            im1 = uint8(arrayfun(@(x) find(abs(d(:)-x) == min(abs(d(:)-x))),current_im));

            % K means
            pixel_labels = imsegkmeans(im1,num_classes);
            
            
            %% FIGURE THIS OUT 
            for n = 2:num_classes


                % Mask the pixel lables
                mask = pixel_labels == n;
                %im2 = im1 < 5 == 0;
                growth_band =  uint8(im1) .* uint8(mask);
                growth_band = imbinarize(growth_band);

                % Make Get the convex image to measure eccentricity 
                convex_image = regionprops(growth_band, 'ConvexImage');
                % Pad the image with 0's to make the measurement easier
               % convex_image.ConvexImage = padarray(convex_image.ConvexImage, [5,5]);
               
               
                if length(convex_image) ~=1
                    continue
                end
                
                % Check to make sure the correct number of peaks were found 
                if length(convex_image) == 1
                    
                    % Produce Convex Image and Measure Eccentricity
                    % convex_image_pixel_list = regionprops('table',convex_image.ConvexImage, 'PixelList', 'Eccentricity');

                    % Measure Feret Dimensions
                    
                    [boundary, ~] = bwboundaries(growth_band);
                    coordinates = [boundary{1}(:,2), boundary{1}(:,1)];
                    % Rotate to orient along major axis
                    rotated = pca_rotate(coordinates);
                    minor_length_temp = abs(min(rotated(:,1))) +  max(rotated(:,1));
                    major_length_temp = abs(min(rotated(:,2))) +  max(rotated(:,2));

                    
                    %feret_table = feretProperties(convex_image_pixel_list);
                    
                    if major_length_temp > 110
                        % Sometimes it measures the background instead of
                        % the ooid :/
                        % If it does this, take the complementary image and
                        % rerun the measurement
                        growth_band = imcomplement(growth_band);
                        
                        [boundary, ~] = bwboundaries(growth_band);
                        coordinates = [boundary{1}(:,2), boundary{1}(:,1)];
                        % Rotate to orient along major axis
                        rotated = pca_rotate(coordinates);
                        minor_length_temp = abs(min(rotated(:,1))) +  max(rotated(:,1));
                        major_length_temp = abs(min(rotated(:,2))) +  max(rotated(:,2));                       
                        
                        if major_length_temp > 110 
                            count1 = count1 +1;
                            keyboard
                            disp(['oops 1 # Bolton is a dumb dumb' num2str(count1)])
                        else
                            % Record Values
                            major_length_t(1, n-1) = major_length_temp;
                            minor_length_t(1, n-1) = minor_length_temp;

                            major_length(j, rep) = max(major_length_t);

                            minor_length(j, rep) = max(minor_length_t);
                        end
                      
                    else

                        % Record Values
                        major_length_t(1, n-1) = major_length_temp;
                        minor_length_t(1, n-1) = minor_length_temp;

                        major_length(j, rep) = max(major_length_t);

                        minor_length(j, rep) = max(minor_length_t);
                        
                       if major_length == 0
                           keyboard
                       end

                      
                    end
                    

                else
%                     keyboard
                    % Sometimes the wrong number of peaks are found 
                    count1 = count1 +1; 
                    disp(['oops #' num2str(count1)])
                end


            end
            
          
        end
        
    else
        % There is only one layer, so the K-Means isn't necessary
%{ 
        for j=1:number_obs

            % Make Get the convex image to measure eccentricity 
            convex_image = regionprops(sliced_imgs{j}, 'ConvexImage', 'Eccentricity');
    
            % Produce Convex Image and Measure Eccentricity
            convex_image_pixel_list = regionprops('table',convex_image.ConvexImage, 'PixelList');


            try
                % Measure Feret Dimensions
                feret_table = feretProperties(convex_image_pixel_list);
                            % Record Values
                major_length(j, rep) = feret_table.MaxFeretDiameter;
                minor_length(j, rep) = feret_table.MinFeretDiameter;
                %eccentricity(j, rep) = convex_image_pixel_list.Eccentricity;
                
            catch ME
                disp('Pesky Error was skipped!')
                continue
            end
        
        end
%}        

    end
    
    
end


