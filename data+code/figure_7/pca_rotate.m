function rotated = pca_rotate(coordinates)
% This function rotates a 2D or 3D data set so that it is oriented by its
% pricipal components
%% IN 
% coordinates:  n x 3 set of datapoints
%% OUT
% rotated: n x 3 set of datapoints after rotation
%
% Example
% coordinates = [0, 0, 0;
%                1, 1, 1;
%                5, 4, 3];
% rotated = pca_rotate(coordinates)
%% Bolton
% April 2019
%% ==============================  Begin  ==============================
if size(coordinates,2) == 2
    x = coordinates(:,1); y = coordinates(:,2);
    
    % Remove the mean
    x0 = x-mean(x); y0 = y-mean(y); 

    % Center the data near the origin 
    P1 = [x0 y0]; 
    
    % Take the pca
    coeff=pca(coordinates);

    % The 3 principal vectors of the input coordinates

    e1=coeff(:,1)'; e2=coeff(:,2)';  

    % 3 unit vectors
    n1=[1 0 ]; n2=[0 1 ];   % 2 unit vector Ox,Oy,Oz

    % Rotation Matrix
    R=[e2;e1];  % rotation matrix 
    
    % Finding the new rotated coordinates
    rotated = (R*P1')'; %new data scaled back to the size of the input coordinates 
    rotated_coordinates = [rotated(:,1)+mean(x),rotated(:,2)+mean(y)];
    
    
    
    
elseif size(coordinates,2) == 3
    % Get coordinates into 3 different vectors
    x = coordinates(:,1) ; y = coordinates(:,2) ;z = coordinates(:,3) ; % get (x,y,z) coordinate  
    
    % Remove the mean
    x0 = x-mean(x) ; y0 = y-mean(y) ; z0 = z-mean(z) ; 

    % Center the data near the origin 
    P1 = [x0 y0 z0]; 

    % Take the pca
    coeff=pca(coordinates);

    % The 3 principal vectors of the input coordinates
    e1=coeff(:,1)'; e2=coeff(:,2)' ;e3=coeff(:,3)';  

    % 3 unit vectors
    n1=[1 0 0]  ; n2=[0 1 0]   ; n3=[0 0 1];   % 3 unit vector Ox,Oy,Oz
    
    % Rotation Matrix
    R=[e2;e1;e3];  % rotation matrix 
    
    % Finding the new rotated coordinates
    rotated = (R*P1')'; %new data scaled back to the size of the input coordinates 
    rotated_coordinates = [rotated(:,1)+mean(x),rotated(:,2)+mean(y),rotated(:,3)+mean(z)];
    
    
    
else
    disp('Check the dimensions of the input. Must be 2D or 3D')
    RETURN
end



% All the code below is for cute graphs
%hold on; scatter3(rotated_coordinates(:,1),rotated_coordinates(:,2),rotated_coordinates(:,3),'r.');
%%Plot the original & rotation 3D data
%figure;scatter3(uniq_boundaries(:,1),uniq_boundaries(:,2),uniq_boundaries(:,3),'b.');
%hold on;scatter3(newdata(:,1),newdata(:,2),newdata(:,3),'r.');
%HA=[min(uniq_boundaries(:,1)) min(uniq_boundaries(:,2)) max(uniq_boundaries(:,3))+1];%just for better visualaztion
%hold on;scatter3(HA(:,1),HA(:,2),HA(:,3),'g.');%just for better visualaztion
%legend('original data(P)','rotated data(newdata)');title({'Plot the rotation matrix 3D point data';'(FINAL RESULT)'});