%% This script plots the equlibrium lines for precipitation and abrasion 
% rates as an ooid grows

% Size
radius = 0:.01:5;
diameter = 2*radius;

% Volume 
volume = (4/3) * pi * radius.^3;
volumetric_rate = volume(2:end) - volume(1:end-1);

% Surface Area 
surface_area = 4*pi*radius.^2;
surface_area_rate = surface_area(2:end) - surface_area(1:end-1);
increased_precip = surface_area + 70;
increased_rate = increased_precip(2:end) - increased_precip(1:end-1);

% Plot
figure('Renderer', 'Painters');
hold on; grid on; box on; axis on; 

% Text sizes
tick_label = 12;
ax = gca;
ax.FontSize = tick_label;
ax.GridAlpha = 0.15;
label_size = 16;
title_size = 20;

color = colormap(lines);
plot(diameter, volume, 'Color', color(1,:))
plot(diameter, increased_precip, 'Color', color(2,:), 'LineStyle', '--')
plot(diameter, surface_area, 'Color', color(2,:), 'LineStyle', '--')


xlabel('Diameter', 'FontSize', label_size)
ylabel('Volumetric Rate', 'FontSize', label_size)
title('Ooid Equilibrium Size', 'FontSize', label_size+2)
legend('Abrasion Rate','Precipitation Rate','Increased Precipitation Rate')

%print(gcf, '-depsc', '-painters', 'precip_vs_abrasion')