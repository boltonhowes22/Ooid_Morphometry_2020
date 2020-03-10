load('growth_history.mat')

% The second ooid has a really small nucleus, which interferes with seeing
% developing confidence intervals, so in this I removed the innermost band
growth_history{2} = growth_history{2}(2:end, :);
%% Plot up the growth histories from Triassic Ooids

a_axis = 1;
b_axis = 2;
c_axis = 3;

for this_ooid = 1:length(growth_history)
    for this_band = 1:size(growth_history{this_ooid}, 1)
        ac{this_ooid}(1, this_band)= growth_history{this_ooid}(this_band,a_axis)/growth_history{this_ooid}(this_band,c_axis);
        bc{this_ooid}(1, this_band) = growth_history{this_ooid}(this_band,b_axis)/growth_history{this_ooid}(this_band,c_axis);
    end
end

%%
g3{1} = growth_history{3};
number_obs = 7500;
[aspect_ratio, major_length, confInt, full_ooid] = triassic_monte_carlo(g3, number_obs);
%%
r = [72 45 57]./255;
g = [29 111 185]./255;
b = [111 142 118]./255;
%microns = 11.2;
[~, example1] = maxk(aspect_ratio{1}(:,1), 2);
[~, example2] = mink(aspect_ratio{1}(:,1), 2);

figure('Renderer', 'Painters');
subplot(1,2,1)
for this_ooid = 1:size(ac, 2)
    growthTrajectories()
    box on; hold on; grid off;
    scatter(ac{this_ooid}, bc{this_ooid}, growth_history{this_ooid}(:,3)/20, 'MarkerEdgeColor',[r(this_ooid) g(this_ooid) b(this_ooid)])
    scatter(ac{this_ooid}, bc{this_ooid}, 'MarkerFaceColor', [r(this_ooid) g(this_ooid) b(this_ooid)], 'MarkerEdgeColor', 'none')
    plot(ac{this_ooid}, bc{this_ooid}, 'Color', [r(this_ooid) g(this_ooid) b(this_ooid)])
   
end
line([0 1], [0 1], 'Color', 'k')

xlim([.3 1])
ylim([.65 1])
   
xlabel('a/c')
ylabel('b/c')


scale = [500 1500 2500 3500 4500];
scalex = [.22 .232 .247 .265 .285];
scaley = [.95 .95 .95 .95 .95]; 
scatter(scalex, scaley,scale/100, 'k')





% Guide Lines: How aspect ratio changes as a function of growth. 
% aspect ratios

ar = [.5, .6, .7, .8, .9, .95];
% Initial major  and minor axis
major_init = repmat(median(major_length{1}(:, 1))-300,  1,length(ar));
minor_init = major_init .* ar;

% Growth 
grow = cumsum(repmat(1, 10000, length(ar))) ;

% Grow the axes
major_grow = major_init + grow;
minor_grow = minor_init + grow;

% Calculate aspect ratio
ar_grow =  minor_grow ./ major_grow ;
%
pbaspect([1 1 1])

subplot(1,2,2)
box on; hold on

plot(major_grow, ar_grow, 'k--')
%plot(growth_history{3}(:, 3), aspect_ratio{3}(example1, :))
%plot(growth_history{3}(:, 3), aspect_ratio{3}(example2, :))

for x = 1:length(major_length{1})
    p = plot(major_length{1}(x, :), aspect_ratio{1}(x,:), 'k');
    p.Color = [r(3) g(3) b(3) 1];
end

plot(major_length{1}(example1(2), :), aspect_ratio{1}(example1(2),:), 'k')
plot(major_length{1}(example2(2), :), aspect_ratio{1}(example2(2),:), 'b')

xlabel('c (microns)')
ylabel('Aspect Ratio')
%xlim([min(growth_history{3}(:,3)), max(growth_history{3}(:,3))])
xlim([median(major_length{1}(:, 1))-300, max(major_length{1}(:))])
ylim([.5 1])
pbaspect([1 1 1])

print(gcf, '-depsc', '-painters', 'triassic_growth_histories')

figure('Renderer', 'Painters');
subplot(1,2,1)
imagesc(full_ooid{1}(:,:,example1(2)))
myColorMap = viridis;
myColorMap(1,:) = 1;
colormap(myColorMap);
axis image 

subplot(1,2,2)
imagesc(full_ooid{1}(:,:,example2(2)))
myColorMap = viridis;
myColorMap(1,:) = 1;
colormap(myColorMap);
axis image
print(gcf, '-depsc', '-painters', 'triassic_traced')

%% Calculate the change in sphericity for all the apparent growth paths. 
% must assume they are prolate
% calculate intermediate axis lenght from the major lenght and aspect ratio
% returned from the monte carlo
intermediate_length{1} = aspect_ratio{1}.*major_length{1};

% calculate surface area and volume
volume{1} = (4/3)*pi*major_length{1}.*intermediate_length{1}.*intermediate_length{1};
surface_area{1} = 4*pi*(((major_length{1} .* intermediate_length{1}).^1.6 + (major_length{1} .* intermediate_length{1}).^1.6 + (intermediate_length{1} .* intermediate_length{1}).^1.6)/3).^(1/1.6);

% calculate sphericity and change in sphericity from first growth band to
% the last growth band. 
sphericity{1} = (pi^(1/3).*(6.*volume{1}).^(2/3))./surface_area{1};
delta_sphericity{1} = sphericity{1}(:,end) - sphericity{1}(:,1);

longb = g3{1}(1,3);
interb = g3{1}(1,2);
shortb = g3{1}(1,1);

longe = g3{1}(end,3);
intere = g3{1}(end,2);
shorte = g3{1}(end,1);

begin_volume_true = (4/3)*pi*longb.*interb.*shortb;
end_volume_true = (4/3)*pi*longe.*intere.*shorte;

begin_surface_area = 4*pi*(((longb .* interb).^1.6 + (longb .* shortb).^1.6 + (interb .* shortb).^1.6)/3).^(1/1.6);
end_surface_area = 4*pi*(((longe .* intere).^1.6 + (longe .* shorte).^1.6 + (intere .* shorte).^1.6)/3).^(1/1.6);

begin_sphericity = (pi^(1/3) * (6*begin_volume_true)^(2/3))/begin_surface_area;
end_sphericity = (pi^(1/3) * (6*end_volume_true)^(2/3))/end_surface_area;
true_delta_sphericity = end_sphericity - begin_sphericity;

figure('Renderer', 'Painters');
box on; grid on; hold on;
histogram(delta_sphericity{1}, 'Normalization', 'Probability')
line([true_delta_sphericity true_delta_sphericity],[0 1])
ylim([0 .14])
print(gcf, '-depsc', '-painters', 'triassic_sphericity_hist')
