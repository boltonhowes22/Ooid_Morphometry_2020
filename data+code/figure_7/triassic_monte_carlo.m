function [aspect_ratio,major_length_rescale, confInt, full_ooid] = growthHistoryConfInt2D(growth_history, number_obs)
% This code is to show the range of growth histories you would get for get
% for the Triassic giant ooids if you only had 2D slices. 
% INNPUT
%{
number_obs: number of slices through the 2D
growth_history{1} = [50 75 80;
                     70 100 105;
                     100 115 120;
                     155 190 200;
                     230 250 270];  

%}      
% OUTPUT
% confInt: 95% conf intervals on growth paths from 2D
% full_ooid: all examples where every growth band is in the slice
%
% Bolton 
% Feb 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% long axis length
maxis_length = 100;

% Rescale the axes to have same proportions, but a maximum axis length of
% 'maxis_length'. Need to rescale because if the meshgrid is too big, it
% will take forever to run
for this_ooid = 1:length(growth_history)
    
    % Rescale
    growth_history_scaled{this_ooid} = maxis_length/max(growth_history{this_ooid}(:)) * growth_history{this_ooid};

end

                 
% Create a meshgrid
A = -maxis_length-5:maxis_length+5;
[X,Y,Z] = meshgrid(A,A,A);

%% OK, let's make the ooids

for this_ooid = 1:length(growth_history_scaled)
    
    % Create a meshgrid for the ooids
    ooid_bands{this_ooid} = zeros(size(A, 2),size(A, 2),size(A, 2));
    
    for this_band = 1:length(growth_history_scaled{this_ooid})
        
        % Is it inside the growth band? 
        band = (X./growth_history_scaled{this_ooid}(this_band, 3)).^2 + (Y./growth_history_scaled{this_ooid}(this_band, 2)).^2 + (Z./growth_history_scaled{this_ooid}(this_band, 1)).^2 <= 1;
        
        % Add it to the existing ooid 
        ooid_bands{this_ooid} = ooid_bands{this_ooid} + band;
        
    end
end

ooid_bands = cellfun(@(x) x*30, ooid_bands, 'un', 0);

%% Now let's slice them 


for this_ooid = 1:length(ooid_bands)
    sliced_imgs(this_ooid, :) = shape_slicer(ooid_bands{this_ooid}, number_obs);
end

%%

num_layers = zeros(number_obs, 1);
if number_obs > 10
    h = waitbar(0,'Working Hard Identifying Number of Layers...');
end


for this_ooid = 1:size(sliced_imgs, 1)
    tracker_complete_layers = 0;
    for this_slice = 1:size(sliced_imgs,2)

        num_all_layers = size(growth_history_scaled{this_ooid},1);
        current_im = sliced_imgs{this_ooid, this_slice};

        % Histogram of pixel values
        p = imhist(uint8(current_im));

        % Find peaks on the histogram
        [pks, ~] = findpeaks(p, 'MinPeakHeight', 100);

        % Number of classes in the kmeans is the number of peaks in the histogram + 1 for the
        % background 
        num_classes = length(pks) + 1;


        % How many layers are being collected? This will be used to test the
        % likelihood of intersecting the nucleus of an ooid
        num_layers(this_ooid,this_slice) = length(pks);
        

        % Pull all samples that intersect all layers
        if length(pks) ==  num_all_layers

            tracker_complete_layers = 1 + tracker_complete_layers;
            full_ooid{this_ooid}(:,:,tracker_complete_layers) = current_im;
        end


        if number_obs > 10
            waitbar(this_band / length(sliced_imgs))
        end

    end
end
close(h)

%% Now, let's segment the images and measure the dimensions of the ellipses

count = 0;

for this_ooid = 1:length(growth_history)
 
    % Number of growth bands
    num_classes = size(growth_history{this_ooid}, 1);

    
   for this_slice = 1:size(full_ooid{this_ooid}, 3)
       current_im = full_ooid{this_ooid}(:,:, this_slice);
       values = unique(current_im);
       values = values(2:end);
       
       for this_band = 1:length(values)
            mask = current_im == values(this_band);
            
            % convert to uint8 and turn the growth band into a binary
            growth_band =  uint8(current_im) .* uint8(mask);
            growth_band = imbinarize(growth_band);
            
            
            % Find just the convex image
            convex_image = regionprops(growth_band, 'ConvexImage', 'Eccentricity');
    
            % Check to make sure the correct number of peaks were found 
            if length(convex_image) == 1
                
                % Get boundaries
                [boundary, ~] = bwboundaries(convex_image.ConvexImage);
                coordinates = [boundary{1}(:,2), boundary{1}(:,1)];
                
                % Rotate to orient along major axis
                rotated = pca_rotate(coordinates);
                
                % Measure and temporarily store the values
                minor_length_temp = abs(min(rotated(:,1))) +  max(rotated(:,1));
                major_length_temp = abs(min(rotated(:,2))) +  max(rotated(:,2));
                
                
                % Record Values
                major_length{this_ooid}(this_slice, this_band) = major_length_temp;
                minor_length{this_ooid}(this_slice, this_band) = minor_length_temp;
                

            else
                % Sometimes the wrong number of peaks are found because of
                % kmeans
               
                count = count +1; 
                disp(['Wrong Number of Peaks #' num2str(count)])
            end

       end
       
        if number_obs > 10
            waitbar(this_slice / length(full_ooid))
        end
       
   end
end
%% Calculate the aspect ratios 


for this_ooid = 1:length(full_ooid)
    major_length{this_ooid} = sort(major_length{this_ooid}, 2);
    minor_length{this_ooid} = sort(minor_length{this_ooid}, 2);
end
aspect_ratio = cellfun(@(x,y) x./y, minor_length,major_length, 'UniformOutput',false);

% output for major_length, return to real dimensions instead of scaled for
% the meshgrid
for this_ooid = 1:length(full_ooid)
    major_length_rescale{this_ooid} = max(growth_history{this_ooid}(:)) .* major_length{this_ooid}/(maxis_length*2);
end

%% 95% confidence intervals
CINT = @(x,p)prctile(x,abs([0,100]-(100-p)/2)); 

for this_ooid = 1:length(growth_history) 
    confInt{this_ooid} = CINT(aspect_ratio{this_ooid}, 95);
end
end





