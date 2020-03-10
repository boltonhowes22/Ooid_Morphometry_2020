function [true_radii, measured_diameter] = manySpheres(num_spheres, mu, sigma)
% This function turns all points within a sphere in a meshgrid to 1 and 
% outstide of the sphere to 0
% IN 
% radius: radius of the sphere
% num_spheres: number of spheres to run the test on 
% mu: average diameter of the sphere
% sigma: standard deviation of diameter
% OUT
% true_radii: the true radius 
% measured_diameter: radius on th eslice 
% Bolton Howes
% January 2019
% =========================================================================
% Resolution of the Meshgrid

x_res = 1;
y_res = 1;
z_res = 1;



% Preallocate Array
true_radii = zeros(num_spheres,1);
measured_diameter = zeros(num_spheres,1);


for n = 1:num_spheres
    
    % Make radii
    true_radii(n) = normrnd(mu,sigma);
    mesh_size = ceil(true_radii(n));
    
    % Build Meshgrid
    x = -mesh_size:x_res:mesh_size;
    y = -mesh_size:y_res:mesh_size; 
    z = -mesh_size:z_res:mesh_size;
    [X, Y, Z] = meshgrid(x,y,z);
    
    % Find if points in sphere
    in_sphere = sphereMaker(true_radii(n), X, Y, Z);
    
    % Random slice
    slice = randi([1 size(in_sphere,3)], 1, 1);
    view_2d = in_sphere(:,:,slice);                     % Random cut through sphere
    diameter = max(sum(view_2d));                       % diameter of circle in cut
    slice_diameter = (x_res * diameter)/2;                % Scale and make into radius
    
    
    % Store
    measured_diameter(n) = slice_diameter;
    
end

end

function in_sphere = sphereMaker(radius, X, Y, Z)
% This function turns all points within a sphere in a meshgrid to 1 and 
% outstide of the sphere to 0
% IN 
% radius: radius of the sphere
% X: X direction of the meshgrid
% Y: Y direction of the meshgrid
% Z: Z direction of the meshgrid
% OUT
% in_sphere: a meshgrid with the sphere as 1 and nonsphere as 0
% Bolton Howes
% January 2019
% =========================================================================


    in_sphere = (X.^2 + Y.^2 + Z.^2 - radius^2) <= radius;
end



