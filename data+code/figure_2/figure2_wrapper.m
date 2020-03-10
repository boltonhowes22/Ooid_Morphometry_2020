%% This script makes all plots for Figure 2
% **Note that in its entirety, this will take several hours to run**

clear all; close all

%% Figure 2A

% Number and distribution of spheres
num_spheres = 10000; 
mu = 100;
sigma = 10;

tic
% Measure from distribution
[diameter, measured_diameter] = manySpheres(num_spheres, mu, sigma);

% Normalize
n_diameter = diameter./max(diameter);
n_measured_diameter = measured_diameter./max(diameter);

% Plot
bin_width = .02;

figure('Renderer', 'painters')
box on; grid on; hold on;
histogram(n_diameter, 'BinWidth', bin_width, 'Normalization', 'Probability')
histogram(n_measured_diameter, 'BinWidth', bin_width, 'Normalization', 'Probability')
ylabel('Probability')
xlabel('Diameter')
hold off;

print(gcf, '-depsc', '-painters', 'Erros-2D-Spheres')
toc

%% Figure 2b
% This code will demonstrate the errors of measuring a distribution of ellipsoidal grains

mu_major = 50;                  % radius of major axis 
mu_intermediate_aspect = .8;
mu_minor_aspect = .6;
sigma_major = 5;
sigma_aspect = .1;


num_layers = 1;                 % 1 growth band
num_obs = 1;                    % number of observations
reps = 2000;                   % not calculating mean/median, so only 1 rep
nuc_bool = 0;                   % no nucleus
nuc_perc = 1;                   % no nucleus % necessary without a nucleus


for x = 1:reps
    
    true_major(x) = normrnd(mu_major,sigma_major);
    true_intermediate(x) = normrnd(mu_intermediate_aspect,sigma_aspect);
    true_minor(x) = normrnd(mu_minor_aspect,sigma_aspect);
    
    % Make the measurements
    [measured_major_length(x), measured_minor_length(x)] = simple_ellipsoid_dist_error(true_major(x), true_intermediate(x),...
        true_minor(x),...
        num_layers,... 
        num_obs,... 
        nuc_bool,... 
        nuc_perc);

end

% Plot
bin_width = .03;

figure('Renderer', 'painters', 'Position', [10 10 700 1350])
subplot(3,2,2)
box on; grid on; hold on;
%line([aspect_ratio_bc*1, aspect_ratio_bc*1], [0, 1], 'LineWidth', 6, 'color', [.1602 .4259, .5586])
%line([aspect_ratio_ac*1, aspect_ratio_ac*1], [0, 1], 'LineWidth', 6, 'color', [.1602 .4259, .5586])
histogram(measured_major_length./max(2.*true_major), 'BinWidth', bin_width,'Normalization', 'Probability', 'FaceColor', [.1016 .6719 .5039])
histogram((2.*true_major)./max(2.*true_major), 'BinWidth', bin_width,'Normalization', 'Probability', 'FaceColor', [.1602 .4259, .5586])
xlim([-.05 1])
%ylim([0 0.15])
xlabel('Major Axis Lengths')
ylabel('Probability')
hold off
print(gcf, '-depsc', '-painters', 'Ellipsoid-Distribution-Error')

%% Figure 2c-d
% This code will demonstrate the errors in measuring a single ellipsoid in 2D


aspect_ratio_ac = .8;          % Only interested in one aspect ratio
aspect_ratio_bc = .60;
num_layers = 1;                 % 1 growth band
num_obs = 20000;                 % number of observations
reps = 30;                       % not calculating mean/median, so only 1 rep
nuc_bool = 0;                   % no nucleus
nuc_perc = 1;                   % no nucleus % necessary without a nucleus

tic
% Make the measurements
[major_length, minor_length] = simple_ellipsoid_error(aspect_ratio_ac,...
    aspect_ratio_bc,...
    num_layers,... 
    num_obs,... 
    nuc_bool,... 
    nuc_perc);
toc
% Make the Graphs
nbins = ceil(1 + log2(numel(1:num_obs))*2);

% Plot
figure('Renderer', 'painters', 'Position', [10 10 700 1350])
subplot(3,2,3)
box on; grid on; hold on;
line([aspect_ratio_bc*1, aspect_ratio_bc*1], [0, 1], 'LineWidth', 6, 'color', [.1602 .4259, .5586])
line([aspect_ratio_ac*1, aspect_ratio_ac*1], [0, 1], 'LineWidth', 6, 'color', [.1602 .4259, .5586])
histogram(minor_length./max(major_length), nbins,'Normalization', 'Probability', 'FaceColor', [.1016 .6719 .5039])
xlim([0 1])
ylim([0 0.3])
xlabel('Minor Axis Lengths')
ylabel('Probability')
hold off

subplot(3,2,4)
box on; grid on; hold on;
line([1, 1], [0, 1], 'LineWidth', 6, 'color', [.1602 .4259, .5586])
histogram(major_length./max(major_length), nbins, 'Normalization', 'Probability', 'FaceColor', [.1016 .6719 .5039])
xlim([0 1.])
ylim([0 0.2])
xlabel('Major Axis Lengths')
ylabel('Probability')
hold off

print(gcf, '-depsc', '-painters', 'Erros-2D-Ellipsoids')

%% Figure 2e
% Error in measuring ellipsoids with different relative axes lengths

aspect_ratio_ac = [.5 .6 .7 .8 .9 1];
aspect_ratio_bc = [.5 .6 .7 .8 .9 1];
num_layers = 1;                 % 1 growth band
num_obs = 100;                  % number of observations
reps = 30;                      % number of reps
nuc_bool = 0;                   % no nucleus
nuc_perc = 1;                   % no nucleus % necessary without a nucleus

major_length = zeros(length(aspect_ratio_ac), length(aspect_ratio_bc), num_obs);
minor_length = zeros(length(aspect_ratio_ac), length(aspect_ratio_bc), num_obs);

h = waitbar(0,'Working...');
j = 0;
tic
true_long_axis = 100;
for z = 1:reps
    for x = 1:length(aspect_ratio_ac)
        true_intermediate_axis(x) = true_long_axis * aspect_ratio_ac(x);
        true_short_axis(x) = true_long_axis * aspect_ratio_bc(x);
        for y = 1:length(aspect_ratio_bc)
        
            [major_length(x,y,:), minor_length(x,y,:)] = simple_ellipsoid_error(aspect_ratio_ac(x),...
                aspect_ratio_bc(y),...
                num_layers,...
                num_obs,...
                nuc_bool,...
                nuc_perc);

            j = j + 1;
            waitbar(j/(length(aspect_ratio_bc)*length(aspect_ratio_ac)*reps))
        end
        
        med_major_length(:,:,z) = median(major_length, 3);
        med_minor_length(:,:,z) = median(minor_length, 3);
        
        
        
        
    end
    
end
toc
close(h)

% Calculate error
avg_perc_error_intermediate_axis = mean((med_minor_length - true_intermediate_axis)./true_intermediate_axis, 3);
avg_perc_error_long_axis = mean((med_major_length - true_long_axis)/true_long_axis, 3);


% Plot
figure('Renderer', 'painters', 'Position', [10 10 700 1350])
subplot(3,2,5)
box on; grid on; hold on;
for x = 1:length(aspect_ratio_ac)
    for y = 1:length(aspect_ratio_bc)
         scatter(aspect_ratio_ac(x), aspect_ratio_bc(y), 200, 100*abs(avg_perc_error_long_axis(x,y)), 'filled')
    end
end
colormap(viridis)
h = colorbar;
caxis([0 100*max(abs(avg_perc_error_long_axis(:)))])
xlim([.45 1.05])
ylim([.45 1.05])
pbaspect([1 1 1])
xlabel('a/c', 'FontSize', 12)
ylabel('b/c', 'FontSize', 12)
print(gcf, '-depsc', '-painters', 'ac-bc-mean-error')

% Uncomment if you want to view the intermediate axis errors
%{
subplot(3,2,6)
box on; grid on; hold on;
for x = 1:length(aspect_ratio_ac)
    for y = 1:length(aspect_ratio_bc)
         scatter(aspect_ratio_ac(x), aspect_ratio_bc(y), 200, 100*abs(avg_perc_error_intermediate_axis(x,y)), 'filled')
    end
end

colormap(viridis)
h = colorbar;
caxis([0 100*max(abs(avg_perc_error_long_axis(:)))])
xlim([.45 1.05])
ylim([.45 1.05])
pbaspect([1 1 1])
xlabel('a/c', 'FontSize', 12)
ylabel('b/c', 'FontSize', 12)
%}

%% Figure 2f
% Caclulate error for different nucleus sizes and ooid shapes
clear all
aspect_ratio_n = [.5, .6, .7, .8, .9, 1];
num_layers = 2; 
num_obs = 100;
reps = 10;
nuc_bool = 1;
nuc_percent = .5;

% To get the volume % of the nucleus, these are the length of the major
% axis of the nucleus

% T .  he lenghts for 5, 10, 15, 20, 25 percent 
percent_of_volume = [5, 10, 15, 20, 25];
nuc_major = [18.4202, 23.20796 26.5664, 29.2402, 31.49081];
major_length_t = zeros(num_obs, reps, length(aspect_ratio_n));
minor_length_t = zeros(num_obs, reps, length(aspect_ratio_n));

tic
h = waitbar(0,'Working...');
for j = 1:length(nuc_major)
 
    for n = 1:length(aspect_ratio_n)

        [major_length_t(:,:,n), minor_length_t(:,:,n)] = double_ellipsoid_error(aspect_ratio_n(n), num_layers, num_obs, reps, nuc_bool, nuc_major(j));    
        
        waitbar((j*n)/(length(nuc_major)*length(aspect_ratio_n)))
    end
    major_length{j} = major_length_t;
    
    minor_length{j} = minor_length_t;
    
end
close(h)
toc

%Calculate Medians: The true answer is 100
 for perc = 1:length(nuc_major)
    a = major_length{perc};
    
    for ar = 1:length(aspect_ratio_n)
        for rep = 1:size(a,2)

            idx = a(:,rep,1) ~= 0;
            m(1,rep, ar, perc) = median(a(idx,rep,ar));

        end
    end
 end
 
% Calculate Error
true = 100;
percent_error = squeeze(abs((true - m))/true * 100);

size(percent_error)
mean_error = squeeze(mean(percent_error,1));

%ß Plot
figure('Renderer', 'painters', 'Position', [10 10 700 1350])
subplot(3,2,6)
percent_of_volume = [5, 10, 15, 20, 25];
oblateness = (true - (aspect_ratio_n.*true))./true;

mean_error_plot = flipud(mean_error);
box on; grid on; hold on;
%c = viridis(ceil(max(mean_error(:))));
for n = 1:5%length(percent_of_volume)
    for j =1:6%length(aspect_ratio_n)
        scatter(oblateness(j), percent_of_volume(n),200, mean_error(j,n), 'filled')
    end
end
colormap(viridis)
colorbar
caxis([0 max(mean_error(:))])
xlabel('Oblateness')
ylabel('Nucleus Percent of Total Volume')
xlim([-.05 .55]) 
ylim([2.5, 27.5])
print(gcf, '-depsc', '-painters', 'nucleus-mean-error')
