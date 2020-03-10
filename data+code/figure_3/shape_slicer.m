function slices = shape_slicer(in_shape, number_obs)
    %% Make a meshgrid
    % Define step size. These can be adjusted to the resolution of the
    % images. 
    x_step = 1;
    y_step = 1;
    z_step = 1;
    
    % Get the size of the input meshgrid
    
    [L, M, S] = size(in_shape);
    
    L = (L-1)/2;
    M = (M-1)/2;
    S = (S-1)/2;
    
    min_L = (-1*L);
    L = min_L:L;
    min_M = (-1*M);
    M = min_M:M;
    min_S = (-1*S);
    S = min_S:S;
    
    [X,Y,Z] = meshgrid(L,M,S);
    
    
    x_solid = X(in_shape ~= 0);
    y_solid = Y(in_shape ~= 0);
    z_solid = Z(in_shape ~= 0);
    

    %% RANDOM SLICING
    scale_factor = 1;

    % Pick a number of points that will serve as the origin of the slice
    point_idxs = datasample(1:length(x_solid), number_obs);
    origin_points = [x_solid(point_idxs), y_solid(point_idxs), z_solid(point_idxs)];

    % Next, we're going to gernerate random numbers on a sphere -> to be
    % treated as normal vectors for a slicing plane
    normal_vectors = produce_random_normals(number_obs);

    % Let's produce a sectional slice that is 2 * max_dim by 2 *
    % second_max_dim
    % First, get the size of hollow volume, and sort descending 
    
    % CHANGES TO AKSHAYS CODE HERE
    size_solid_volume = sort([max(x_solid(:)) - min(x_solid(:)), ...
        max(y_solid(:)) - min(y_solid(:)), ...
        max(z_solid(:)) - min(z_solid(:))], 'descend');
    [slice_x, slice_y] = meshgrid(-size_solid_volume(1):scale_factor:size_solid_volume(1), ...
        -size_solid_volume(1):scale_factor:size_solid_volume(1));
    
    % Zeros z
    slice_z = zeros(size(slice_x));

    % Four corners of surface (for slice visualization)
    min_slice_x = min(slice_x(:));
    max_slice_x = max(slice_x(:));
    min_slice_y = min(slice_y(:));
    max_slice_y = max(slice_y(:));
    four_corners = [min_slice_x, min_slice_y, 0; ...
        min_slice_x, max_slice_y, 0; ...
        max_slice_x, max_slice_y, 0; ... 
        max_slice_x, min_slice_y, 0; ...
        ];
    % Declare normal vector
    slice_normal = [0 0 1];

    % Gridded interpolant
    p = [2 1 3];
    x_p = double(permute(X,p));
    y_p = double(permute(Y,p));
    z_p = double(permute(Z,p));
    vol_p = double(permute(in_shape, p));
    gridded_func = griddedInterpolant(x_p, y_p, z_p, vol_p,  'nearest');
    diff_collector = 0;
%     
    if number_obs > 10
        %h = waitbar(0,'Working Hard Slicing and Dicing...');
    end
    
    for pt = 1:number_obs

        % First, the origin and normal vector of this pt
        this_origin = origin_points(pt, :);
        this_normal = normal_vectors(pt, :);

        % Next, figure out rotation from slice_normal to this_normal
        rot = vrrotvec(slice_normal, this_normal);

        % Define a surface
        %slice_surf = surf(slice_x, slice_y, slice_z);
        slice_surf = surf(slice_x + this_origin(1), slice_y  + this_origin(2), slice_z + this_origin(3));

        % Rotate this surface
        rotate(slice_surf, rot(1:3), rad2deg(rot(4)), this_origin); 
        slice_surf_x = slice_surf.XData;
        slice_surf_y = slice_surf.YData;
        slice_surf_z = slice_surf.ZData;
        delete(slice_surf);
        
        sliced_interp = gridded_func(slice_surf_x, slice_surf_y, slice_surf_z);
        sliced_interp_surf = surf(slice_surf_x, slice_surf_y, sliced_interp);
        sliced_interp_img = sliced_interp_surf.CData;
        delete(sliced_interp_surf);

        sliced_interp_img(isnan(sliced_interp_img)) = 0;
        % Closing holes
        structuring = strel('disk',5);

        % Close it
        sliced_interp_img = imclose(sliced_interp_img, structuring);
        
        slices{pt} = sliced_interp_img;
        
        if number_obs > 10
            %waitbar(pt / number_obs)
        end
        
        % slice_collector(:,:,pt) = sliced_interp_img;
        %{
        % Visualizing sections
        axis image
        % Convert rot rotation into a matrix
        rot_m = vrrotvec2mat(rot);
        % Apply matrix rotation to four corners
        rot_corners = four_corners * rot_m';
        % Translate
        rot_corners = rot_corners + repmat(this_origin, 4, 1);
        % Add first row to bottom (in order to close the plane
        rot_corners = [rot_corners; rot_corners(1,:)];
        plot3(rot_corners(:,1), rot_corners(:,2), rot_corners(:,3));
        patch(rot_corners(:,1), rot_corners(:,2), rot_corners(:,3), [0 0 0], 'FaceAlpha',.25);
        grid on 
        %}
    end 
        %imshow(sliced_interp_img)
        
%         %%
%         p = imhist(sliced_interp_img);
%         [pks, ~] = findpeaks(p, 'MinPeakHeight', 200);
%         num_classes = length(pks);
%         
%         d=linspace(min(sliced_interp_img(:)),max(sliced_interp_img(:)),256);
%         im1=uint8(arrayfun(@(x) find(abs(d(:)-x)==min(abs(d(:)-x))),sliced_interp_img));
%         pixel_labels = imsegkmeans(im1,num_classes);
%         % B = labeloverlay(sliced_interp_img, pixel_labels); 
%         
%         
%     
%         
%         mask3 = pixel_labels == 2;
%         im2 = im1 < 30 == 0;
%         cluster3 =  (uint8(im2)*255) .* uint8(mask3);
%         cluster3 = im2bw(cluster3);
%         imshow(cluster3)
%         test = regionprops(cluster3, 'FilledImage', 'ConvexImage', 'MajorAxisLength', 'MinorAxisLength', 'Solidity', 'Image')
%         
%         
%         % Write a For loop that fills each ring created by the mask
%         % Then calculate eccentricity and the aspect ratio. 
%         
%         
%         % Stats
%         stat = regionprops(sliced_interp_img, 'MajorAxisLength', 'MinorAxisLength', 'Area', 'Eccentricity');
%         area(pt) = stat.Area;
%         major_axis_length(pt) = stat.MajorAxisLength;
%         minor_axis_length(pt) = stat.MinorAxisLength;
%         
%         eccentricity(pt) = stat.Eccentricity;
        
        % Save Slices
        % save(fullfile(storage, ['slice_', num2str(pt), '.mat']), 'sliced_interp_img');
    if number_obs > 10
        %close(h)
    else
    end
end



%%
function [normal_vectors] = produce_random_normals(num)
    % Taken from https://www.mathworks.com/help/matlab/math/numbers-placed-randomly-within-volume-of-sphere.html
    % Calculate elevation angles -> NOT UNIFORM
    rvals = 2*rand(num,1)-1;
    elevation = asin(rvals);
    % Calculate azimuth angles -> UNIFORM
    azimuth = 2*pi*rand(num,1);
    % Radius of 1
    radius = 1;
    [x,y,z] = sph2cart(azimuth,elevation,radius);
    normal_vectors = [x,y,z];
end