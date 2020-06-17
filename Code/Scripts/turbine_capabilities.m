%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:

%%This script is used to plot the capabilities of the different turbines.
%%This includes the capability curves and the active/reactive power outputs
%%at different windspeeds. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%%Generate wind speeds
stepsize = 0.1;
vmin = 0;
vmax = 30;
v = vmin:stepsize:vmax;

%%Compute P and Q
[P,Q] = compute_pq_wtg(v);

%%Stuff for plot
endsamp_string = 25/stepsize + 1;

for i = 1:length(P(1,:))
figure(1)
subplot(4,4,i)
plot(P(1:endsamp_string,i),Q(1:endsamp_string,i))    
xlabel('P [MW]')
ylabel('Q [MVAr]')

figure(2)
subplot(4,4,i)
plot(v,P(:,i))
xlabel('windspeed [m/s]')
ylabel('P [MW]')

figure(3)
subplot(4,4,i)
plot(v,Q(:,i))
xlabel('windspeed [m/s]')
ylabel('Q [MVAr]')
end