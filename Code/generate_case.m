%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to generate the reactive power boundaries for the
%%WTG and PVG strings as a function of the windspeed and the irradiance. 
%%Furthermore it updates the casefile with the corresponding active power
%%outputs at the strings.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Q_wt_min, Q_wt_max, Q_pv_min, Q_pv_max] = generate_case(windspeed,irradiance)
global CONSTANTS Systemdata Optimisation

%%Compute P and Q of the wtg strings
[P_wtg, Q_wtg] = compute_pq_wtg(windspeed);

%%Symmetric boundaries
Q_wt_max = Q_wtg;
Q_wt_min = -Q_wt_max;

%%Updates the active power outputs of the wtg strings
Systemdata.mpc.gen(Systemdata.wtg_pos,CONSTANTS.PG) = P_wtg;

%%Repeat for pv if the input is existent
if nargin > 1
    %%compute P and Q of the pvg strings
    [P_pvg, Q_pvg] = compute_pq_pvg(irradiance,Optimisation.Npv);
    
    %%Only positive Q
    Q_pv_max = Q_pvg;
    Q_pv_min = zeros(1,length(Q_pvg));

    %%Update casefile
    Systemdata.mpc.gen(Systemdata.pvg_pos,CONSTANTS.PG) = P_pvg;
    
    %%Update bus status for the PV strings
    pv_busvec = [11,16,22,27];
    Systemdata.mpc.bus(pv_busvec,CONSTANTS.BUS_TYPE) = 1;
    
    %%Update branch status for the PV strings
    pv_branchvec = [10,16,23,29];
    Systemdata.mpc.branch(pv_branchvec,CONSTANTS.BR_STATUS)= 1;
else
    Q_pv_min = 0;
    Q_pv_max = 0;
end
end