[trower2018, text1, ~] = xlsread('Trower2018_dataSupplement.xlsx');
%% I had to impute vaules for the sorting of the Shoal Non Crest 
for n = 2:size(text1,1)
    sample_names{n-1} = text1{n, 1};
    environments{n-1} = text1{n, end};
end

for n = 1:size(text1,2)
    measured{n} = text1{1,n};
end
%%
% Isolate non location/Distance measurements
data = trower2018(:, 1:6);
metrics = measured(2:8);

% Numeric identifier of the environment -- makes it easier to iterate
% through environments while plotting
env_num = trower2018(:,end);
%% normalize the data by subtracting the mean and dividing by the standard deviation
for n = 1:size(data,2)
    norm_data(:,n) = (data(:,n) - mean(data(:,n)))/std(data(:,n));
end

%% PCA
[coeff,score,latent,tsquared,explained] = pca(norm_data);
projected_data = norm_data * coeff;

%% V are the eigenvectors, and...
[U,S,V] = svd(norm_data);

for i = 1:size(norm_data(:,1:5),2)
pc{i} = sprintf('%+3.1fD10%+3.1fD50%+3.1fD90%+3.1fSort%+3.1fSpher%+3.1fAsR%+3.1fRound',...
    V(:,i));
end
%% Plot PCA

% Numeric identifier of the environment -- makes it easier to iterate
% through environments while plotting
env_num = trower2018(:,end);

%%
figure('Renderer', 'Painters');
subplot(1,2,1)
environment_array = [{'Foreshore'}, {'Shoal Crest'}, {'Shoal Non-Crest'}, {'Algae-Stabilized'}, {'Inter-Cay Shoal'}];
unique_env = unique(env_num);
hold on
grid on
c = colormap(lines);
c_index = [1,2,2,4,5];
box on 

for n = 1:length(unique_env)
    
    % Combined the Shoal Crest and 2017 Shoal Crest from Trower 2018
    scatter(projected_data(env_num == unique_env(n),1),...
        projected_data(env_num == unique_env(n),2),...
         'SizeData', 130, 'MarkerFaceColor', c(c_index(n),:),...
        'MarkerEdgeColor', 'none')
    
    text(2, -.6 - (.4*n), sprintf('%s', environment_array{n}), 'Color', c(c_index(n),:), 'Fontsize', 14, 'Interpreter', 'Latex')
        
end
tick_label = 12;
ax = gca;
ax.FontSize = tick_label;
ax.GridAlpha = 0.15;

label_size = 16;
title_size = 20;
xlabel(['PC1 Score' sprintf('%3.0f%% ',explained(1)) 'of variance'], 'Interpreter', 'Latex')
ylabel(['PC2 Score' sprintf('%3.0f%% ',explained(2)) 'of variance'], 'Interpreter', 'Latex')


subplot(1,2,2)
h = biplot(coeff(:,1:2),'varlabels',{'D10','D50','D90','Sorting','Sphericity','Aspect Ratio'});
box on
%%
figure('Renderer', 'Painters');
subplot(2,2,2)
bar(explained, 'FaceColor', c(c_index(1),:))
hold on; grid on; box on
ax = gca;
ax.FontSize = tick_label;
ax.GridAlpha = 0.15;

xlabel('Principal Component')
ylabel('% Variance Explained')


subplot(2,2,3)
scatter(V(1,:), V(2,:), 'filled')
ax = gca;
ax.FontSize = tick_label;
ax.GridAlpha = 0.15;

% Axis through Center
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on

xlabel('PC1 Loadings')
ylabel('PC2 Loadings')

subplot(2,2,4)
scatter(V(2,:), V(3,:), 'filled')

% Axis through Center
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
grid on

xlabel('PC2 Loadings')
ylabel('PC3 Loadings')
