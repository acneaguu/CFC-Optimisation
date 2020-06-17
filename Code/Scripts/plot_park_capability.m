%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This script is used to plot the reactive power capability of the WF as
%%a function of the wind speed. It also includes various considered
%%setpoints.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%%Load powerflow data of system
load('workspace_powerflow.mat','results_inject_t');
load('workspace_powerflow.mat','results_loss_im');

%%Calculate the average total branch injection
mean_samples = mean(results_inject_t,2);
tot_inj = sum(mean_samples);

%%Generate wind speeds
vmin = 0;
vmax = 25;
stepsize = 0.1;
v = vmin:stepsize:vmax;

%%Compute the maximum P/Q production at each V
for i = 1:length(v)
    [P,Q] = compute_pq_wtg(v(i));
    capability(i) = sum(-1*Q);
    p_out(i)= sum(P);
end

%% Plot the active power capability
figure(1)
axes_fontsize = 15;
titlesize = 20;
plot(v,p_out);
title('Windfarm Active Power Capability','FontSize',titlesize);
ylabel('P [MW]','FontSize',axes_fontsize);
xlabel('Windspeed[m/s]','FontSize',axes_fontsize)

%% Plot reactive power capabiltiy 
figure(2)
axes_fontsize = 15;
titlesize = 20;

xlabel('Windspeed [m/s]','FontSize',axes_fontsize);
ylabel('Q [MVAr]','FontSize',axes_fontsize);
title('Windfarm Reactive Power Capability','FontSize',titlesize)
ax = gca;
ax.FontSize = axes_fontsize;

%%Plot the park capability with reactor and branch injections
hold on;
plot(v,capability-12,'--b');
legendstr{1} = 'Q_{min} with reactor';
plot(v,capability-12+tot_inj,'b');
legendstr{2} = 'Q_{min} with reactor and branch injections';
plot(v,-capability,'--r');
legendstr{3} = 'Q_{max}';
plot(v,-capability+tot_inj,'r');
legendstr{4} = 'Q_{max} with injections';

%%Plot the setpoints
Q = [-100 -50 0 50 100];
threshold = 12.5/2;

for i = 1:length(Q)
    %%Alternating colours for the setpoints
    if mod(i,2)
        colorsetpoint = '#04cc82';
        colortolerance = '#013220';
    else
        colorsetpoint = '#e59400';
        colortolerance = '#b27300';
    end
    
    %%Horizontal lines for the setpoints
    yline(Q(i),'Color',colorsetpoint,'LineWidth',1.5 )
    yline(Q(i)+threshold,'--','Color',colortolerance,'LineWidth',1.5)
    yline(Q(i)-threshold,'--','Color',colortolerance,'LineWidth',1.5')
    
    %%Save the name for the legend
    N = length(legendstr);
    name = sprintf('%d MVAr',Q(i));
    legendstr{N+1} = strcat("Q_{setpoint} = " + '',name);
    legendstr{N+2} = 'Q_{setpoint} + 6.25 MVAr';
    legendstr{N+3} = 'Q_{setpoint} - 6.25 MVAr';
end

%%Add labels and legend to the figure
lgd = legend(legendstr);
lgd.NumColumns = 3;
lgd.FontSize = 10;