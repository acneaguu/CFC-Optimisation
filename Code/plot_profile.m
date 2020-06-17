%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This script is used to plot the test profile. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%%Generate the test profile: for each wind speed and solar irradiance, 5
%%Qsetpoints are requested
Qref.setpoint =  [-0.286; -0.143; 0; 0.143; 0.286]; 
v = [4.5 4.5 4.5 4.5 4.5 5 5 5 5 5 7 7 7 7 7 12 12 12 12 12 15 15 15 15 15]';
irradiance = [0 0 0 0 0 340 340 340 340 340 680 680 680 680 680 ...
        510 510 510 510 510 170 170 170 170 170];

%%Put the data in a matrix cases    
cases(:,1) = v;
cases(:,2) =repmat(Qref.setpoint,5,1);
cases(:,3) = irradiance;

%%Repeat last entry for legibility
cases(26,1) = 15;
cases(26,2) = 0.286;
cases(26,3) = 170;

%%Compute the P/Q for the different cases in p.u.
%%WTG
[P_wtg, Q_wtg] = compute_pq_wtg(cases(:,1));
P = sum(P_wtg,2)/350;
Q = sum(Q_wtg,2)/350;

%%PVG
[Ppvg,Qpvg] = compute_pq_pvg(cases(:,3),4);
Ppvg = sum(Ppvg,2)/350;
Qpvg = sum(Qpvg,2)/350;

%%Plot variables
axes_fontsize = 15;
titlesize = 20;
Ncase = 1:length(cases(:,1));

%%Plot
fig = figure(1);
set(fig,'defaultAxesColorOrder',[0, 0, 0;0, 0, 0]);
% title('Test Profile Optimisation Unit','FontSize',titlesize)
ax = gca;
ax.FontSize = axes_fontsize;
hold on

%%Left plot is P of the WTG and PVG
yyaxis left
plot(Ncase,P,'-b','LineWidth',1.5);
plot(Ncase,Ppvg,'-','Color','#ff9d3b','LineWidth',1.5);

%%Right plot is Qsetpoint
ylabel('P [p.u.]','FontSize',axes_fontsize)
yyaxis right
stairs(Ncase,cases(:,2),'Color','green');
ylabel('Q [p.u.]','FontSize',axes_fontsize)
xlabel('Case','FontSize',axes_fontsize)

%%Add legend
lgd = legend('P_{available} of all the WTG strings',...
    'P_{available} of all the PVG strings','Q_{requested}');
lgd.Location = 'northwest';
lgd.FontSize = axes_fontsize;