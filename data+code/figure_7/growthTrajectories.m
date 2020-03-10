function growthTrajectories()
timesteps = 5000;
initial_ratio_BC = [.1 .2 .3 .4 .5 .6 .7 .8 .9 1];
initial_ratio_AC = [.0 .0 .0 .0 .0 .0 .0 .0 .0 .0];

for x = 1:length(initial_ratio_BC)
    
    [long_axis(:,x) intermediate_axis(:,x) short_axis(:,x)] = grower(initial_ratio_AC(x),...
                                                      initial_ratio_BC(x),...
                                                      timesteps);
                                               
end
Eac = short_axis./long_axis;
Ebc = intermediate_axis./long_axis;
%{
fig_name = 'growth-trajectories';
f1 = figure(1);
f1.Renderer =  'Painters';
f1.Name = fig_name;
clf
hold on; box on; grid on;
%}

for x = 1:length(initial_ratio_BC)
    p = plot(Eac(:,x), Ebc(:,x), '--');
    hold on 
    p.Color = [0.8 0.8 0.8];
end

plot([0 1], [0 1], 'k')





end
%% Function for Growing 3D Ooids
function [long_axis intermediate_axis short_axis] = grower(initial_ratio_AC, initial_ratio_BC, timesteps)
% Grows an ooid in 3D
% C: long axis 
% B: intermediate axis
% A: short axis
%
% IN 
% initial_ratio_CA: intial c-axis/a-axis ratio
% initial_ratio_BA: intial b-axis/a-axis ratio
% timesteps: number of time steps
% 
% OUT 
% final axis lengths
%
% Bolton Howes
% January 2020

long_axis = 100;
intermediate_axis = long_axis * initial_ratio_BC;
short_axis = long_axis * initial_ratio_AC;
 

for x = 1:timesteps
    long_axis(x+1) = long_axis(end) + 1;
    intermediate_axis(x+1) = intermediate_axis(end) + 1;
    short_axis(x+1) = short_axis(end) + 1;
end


end

