%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to update the casefile with a new solution.
%%If mode is 1, then the Q of the strings, the transformer taps and
%%reactors can be updated. If mode is 2, the active power outputs of the
%%strings are updated. The reactive power outputs of the WTG strings are
%%always updated. The reactive power outputs of the PVG strings, the
%%transformer tap positions and the reactor status. Mode 3 is legacy
%%control for system41 topology.

%%NOTE: initialisation of active power is currently done in generate_case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function update_casefile(Xin,mode)
global Systemdata CONSTANTS Optimisation

%%Updates the controllable variables in the case file
if mode == 1
    %%Update the WTG strings
    Systemdata.mpc.gen(Systemdata.wtg_pos,CONSTANTS.QG:CONSTANTS.QMIN) = ...
        repmat(transpose(Xin(Optimisation.wtg_pos)),1,3);
    
    %%Updates the pvg strings
    Systemdata.mpc.gen(Systemdata.pvg_pos,CONSTANTS.QG:CONSTANTS.QMIN) = ....
        repmat(transpose(Xin(Optimisation.pvg_pos)),1,3);
    
    %%Updates transformer tap positions
    Systemdata.mpc.branch(Systemdata.trans,CONSTANTS.TAP) = Xin(Optimisation.tr_pos);
    
    %%Updates reactor status: bustype and branch status
    Systemdata.mpc.bus(Systemdata.shunts,CONSTANTS.BUS_TYPE) = ...
        4-3*Xin(Optimisation.r_pos); 
    Systemdata.mpc.branch(Systemdata.shuntbranch,CONSTANTS.BR_STATUS) = ...
         Xin(Optimisation.r_pos);

%%Updates the active power     
elseif mode == 2    
    Systemdata.mpc.gen(index_wtg,[2 9:10]) = repmat(transpose(Xin),1,3);

%%Legacy for system 41    
elseif mode == 3    
    Systemdata.mpc.bus(24:end,4) = Xin(np,18).';
    Systemdata.mpc.branch(1,9) = Xin(np,19); %change tf ratio 
    Systemdata.mpc.branch(13,9) = Xin(np,20);
    Systemdata.mpc.bus(2,CONSTANTS.BS) = Xin(np,21);%Changes inductor
    Systemdata.mpc.bus(5,CONSTANTS.BS) = Xin(np,22);%Changes capacitor
end
end