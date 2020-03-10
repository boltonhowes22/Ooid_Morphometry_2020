%axes_3 = load('all_axes.mat');


% Data for E038_PDN 
E038_PDN = [32,0,0,26,115,126,135,410,986,1450,4064,11001,32035,44104,33695,12040,3098,701,288,147,50,11,3,2,0,1,0];
E038_binleft = [0, 25, 30, 37, 45, 55,68,83,102,125,153,188,230,283,346,425,521,639,783,960,1177,1443,1770,2170,2661,3262,4000];
E038_binright = [25, 30, 37, 45, 55,68,83,102,125,153,188,230,283,346,425,521,639,783,960,1177,1443,1770,2170,2661,3262,4000, 5000];
E038_bincenters = (E038_binleft + E038_binright)./2;
%%
C = axes(:,1);
B = axes(:,2);
A = axes(:,3);

figure('Renderer', 'Painters');
subplot(2,2,1)
hold on; box on; grid on;
%ouptut_path
heatscatter(A./C, B./C,[], 'heatscatter_plot')
xlim([0 1]);
ylim([0 1]);
line([0 1],[0 1], 'Color', 'k')
grid on; box on; 
colormap(viridis)
print(gcf, '-depsc', '-painters', 'heatscatter_plot')
%%
subplot(2,2,2)
hold on; box on; grid on;
plot(E038_bincenters, E038_PDN./sum(E038_PDN), 'k')
xlim([0 1000])
num_bins = 10;

edges = [E038_binleft E038_binright(end)];
h = histogram(B, 'BinEdges',edges, 'Normalization', 'probability');
p = histcounts(B, 'BinEdges',edges, 'Normalization', 'probability');

plot(E038_bincenters, p, 'r')
xlabel('Intermediate Axis Length')
ylabel('Probability')
print(gcf, '-depsc', '-painters', 'size_distributions')

%%
sorting = [(prctile(B, 90) - prctile(B, 10))/prctile(B, 50)];
avg_shoreface_sorting = 1.458;
