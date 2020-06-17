%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%-------------------------------------------------------------------------
%%THIS VERSION IS USED FOR THE DETERMINATION OF THE WEIGHTS. 
%%IN THIS VERSION THE OBJECTIVES ARE IN EURO's €\
%%-------------------------------------------------------------------------

%%This function is used to compute the value of the OF of the current
%%solution Xin at t considering also past positions of the transformer taps
%%and the reactor connection. This is done by assigning a cost in euro's to
%%the different objectives. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OF = compute_costs_v2(Xin)
global CONSTANTS Systemdata PFresults Optimisation Results; 
    
%% Ploss
%%Assume 8 cents/kwh = 80 euro/MWh

%%Branch losses
[losses] = get_losses(PFresults);
Ploss_branch = sum(real(losses));

%%Shunt losses (in the busses)
Ploss_shunt = sum(PFresults.bus(:,CONSTANTS.VM) .^ 2 .* PFresults.bus(:,CONSTANTS.GS)); %%shunt absorption real

%%Total losses: 
Ploss_tot = (Ploss_branch + Ploss_shunt)*Optimisation.c1 * Optimisation.timeinterval;
%% Tap switches
%Assume 0.40 euro / tapswitch such that 1 switch = 5kWh

%%Only compute tap switch costs when transformer taps are controlled
if Optimisation.Ntr ~= 0
    %%Tap changes in ratios
    tap_changes_ratio = abs(Xin(Optimisation.tr_pos)-Results(Optimisation.t-1).best_run_solution(Optimisation.tr_pos));

    %%Convert the tap changes from rations to tap positions
    tap_changes = sum((tap_changes_ratio./Systemdata.trstep))*Optimisation.c2;
else
    tap_changes = 0;
end
%% Reactor
%%Assume 0.40 euro/reactor switch such that 1 switch = 5 kWh

%%Only compute reactor costs when reactors are controlled
if Optimisation.Nr ~= 0
    %%Relative reactor change
    reactor_changes = Optimisation.c3*sum(abs(Xin(Optimisation.r_pos) - ...
    Results(Optimisation.t-1).best_run_solution(Optimisation.r_pos))); 
else
    reactor_changes = 0;    
end
%% Extremeness of the setpoints
%Assume 0.4 euro/0.7 mean reactive power distance s.t. 0.7 mean
%distance = 5 kWh

%%Compute the extremeness of the setpoints w.r.t. the upper boundaries of
%%the system. 
extremeness_setpoints = sum(abs(Xin(Optimisation.wtg_pos | Optimisation.pvg_pos))...
/Systemdata.ub(Optimisation.wtg_pos | Optimisation.pvg_pos));

%%Take the average to compute the costs.
extremeness_setpoints = (extremeness_setpoints/(Optimisation.Nturbines+Optimisation.Npv))*Optimisation.c4;
%% Calculate OF
OF = Optimisation.w1*Ploss_tot+Optimisation.w2*tap_changes+...
    Optimisation.w3*reactor_changes + Optimisation.w4*extremeness_setpoints;
end