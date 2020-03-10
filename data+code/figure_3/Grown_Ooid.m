%% DIMENSIONS
% 1
long_axis1 = 100;
medium_axis1 = long_axis1 * .8;
short_axis1 = long_axis1  * .8;

% 2
long_axis2 = 80;
medium_axis2 = long_axis2 * .75;
short_axis2 = long_axis2 * .75;


% 3
long_axis3 = 65;
medium_axis3 = long_axis3 * .7;
short_axis3 = long_axis3 * .7;


% 4
long_axis4 = 55;
medium_axis4 = long_axis4 * .65;
short_axis4 = long_axis4 * .65;

%  5
long_axis5 = 30;
medium_axis5 = long_axis5 * .6;
short_axis5 = long_axis5 * .6;


% Create Meshgrid
A = -long_axis1-5:long_axis1+5;
[X,Y,Z] = meshgrid(A,A,A);

%% In the ellipsoid?
in_ellipsoid1 = (X./long_axis1).^2 + (Y./medium_axis1).^2 + (Z./short_axis1).^2 <= 1;
in_ellipsoid2 = (X./long_axis2).^2 + (Y./medium_axis2).^2 + (Z./short_axis2).^2 <= 1;
in_ellipsoid3 = (X./long_axis3).^2 + (Y./medium_axis3).^2 + (Z./short_axis3).^2 <= 1;
in_ellipsoid4 = (X./long_axis4).^2 + (Y./medium_axis4).^2 + (Z./short_axis4).^2 <= 1;
in_ellipsoid5 = (X./long_axis5).^2 + (Y./medium_axis5).^2 + (Z./short_axis5).^2 <= 1;


in_ellipsoid1 = in_ellipsoid1 .* 50;
in_ellipsoid2 = in_ellipsoid2 .* 50;
in_ellipsoid3 = in_ellipsoid3 .* 50;
in_ellipsoid4 = in_ellipsoid4 .* 50;
in_ellipsoid5 = in_ellipsoid5 .* 50;


in_ellipsoid = (in_ellipsoid1 + in_ellipsoid2 + in_ellipsoid3 + in_ellipsoid4 + in_ellipsoid5)./255;
figure(1234)
imshow(in_ellipsoid(:,:,100))
colorMap = (viridis(256));

% Make background white -- ones(1,1) or black --- zeros(1,1) :
colorMap(1,:) = ones(1,1);
colormap(colorMap);

%%
number_obs = 500;
sliced_imgs = shape_slicer(in_ellipsoid, number_obs);

figure(893524)
montage(sliced_imgs(1:49), 'BorderSize', [1 1])
colormap(colorMap);
title('Random Slices', 'FontSize', 24)


%%

num_layers = zeros(number_obs, 1);
if number_obs > 10
    h = waitbar(0,'Working Hard Identifying Number of Layers...');
end

tracker_complete_layers = 0;

for n = 1:length(sliced_imgs)
    current_im = sliced_imgs{n};
    
    % Histogram of pixel values
    p = imhist(current_im);
    
    % Find peaks on the histogram
    [pks, ~] = findpeaks(p, 'MinPeakHeight', 200);
    
    % Number of classes in the kmeans is the number of peaks in the histogram + 1 for the
    % background 
    num_classes = length(pks) + 1;
    
    % How many layers are being collected? This will be used to test the
    % likelihood of intersecting the nucleus of an ooid
    num_layers(n) = length(pks);
    
    
    % Pull all samples that intersect all layers
    if num_layers(n) ==  5
        tracker_complete_layers = 1 + tracker_complete_layers;
        five{tracker_complete_layers} = sliced_imgs{n};
    end


    if number_obs > 10
        waitbar(n / length(sliced_imgs))
    end

end
close(h)
%%
num_classes = 6;
if number_obs > 10
    h = waitbar(0,'K Means...');
end

% preallocate size of storage
major_length = zeros(length(five), num_classes - 1);
minor_length = zeros(length(five), num_classes - 1);
eccentricity = zeros(length(five), num_classes - 1);

% Track mistaken number of peaks
count = 0;

for j = 1:length(five)
    current_im = five{j};
    
    % Convert image to uint8
    d = linspace(min(current_im(:)), max(current_im(:)),256);
    im1 = uint8(arrayfun(@(x) find(abs(d(:)-x) == min(abs(d(:)-x))),current_im));
    
    % K means
    pixel_labels = imsegkmeans(im1,num_classes);
    
    for n = 2:num_classes

        % Get growth band mask and binarize
        mask = pixel_labels == n;
        im2 = im1 < 5 == 0;
        growth_band =  uint8(im1) .* uint8(mask);
        growth_band = imbinarize(growth_band);
        

        convex_image = regionprops(growth_band, 'ConvexImage', 'Eccentricity');
        
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

                % Temporarily record
                minor_length_temp = abs(min(rotated(:,1))) +  max(rotated(:,1));
                major_length_temp = abs(min(rotated(:,2))) +  max(rotated(:,2));
                
            

            
            % Record Values
            major_length(j, n-1) = major_length_temp;
            minor_length(j, n-1) = minor_length_temp;


        else
            % Sometimes the wrong number of peaks are found 
            count = count +1; 
            disp(['oops #' num2str(count)])
        end
        

    end
    
    if number_obs > 10
        waitbar(j / length(five))
    end
end
%% Sort
% Now there are 5 columns, with the metric for each growth band in the
% column 
major_length = sort(major_length,2);
minor_length = sort(minor_length,2);

% Rounding erros in the measuremnt + the linear interpolation occasionally 
% makes a measure slightly over 200 (only off by 1-ish pixel maximum) 
major_length(major_length > 2* long_axis1) = 2* long_axis1;


%% Plot
label_size = 12; 

% Boxplot of (a=b)/c
ac1 = medium_axis1/long_axis1;
ac2 = medium_axis2/long_axis2;
ac3 = medium_axis3/long_axis3;
ac4 = medium_axis4/long_axis4;
ac5 = medium_axis5/long_axis5;


flattening = minor_length./major_length;

figure('Renderer', 'Painters'); hold on;
boxplot(flattening)
scatter(1, ac5, 200, 'k', 'filled')
scatter(2, ac4, 200, 'k', 'filled')
scatter(3, ac3, 200, 'k', 'filled')
scatter(4, ac2, 200, 'k', 'filled')
scatter(5, ac1, 200, 'k', 'filled')
xlabel('Growth Band', 'FontSize', label_size)
ylabel('(a=b)/c', 'FontSize', label_size)

ylim([0.5 1])

grid on
