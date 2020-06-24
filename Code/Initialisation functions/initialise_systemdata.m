%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to initialise the system topology and powerflow
%%options. The desired topology can be entered with the input topology.
%%This function also initialises logic vectors used to change the
%%generation of the WTG and PVG strings. For this, ensure that the
%%following format is present in the generator variable of MATPOWER: 
%%-Firstly, the slack bus.
%%-Secondly, Nwtg WTG strings.
%%-Lastly, Npv PVG strings.
%%Also, each non-zero ratio brach is considered to be a controllable
%%transformer and all busses with a shunt susceptance are considered as a
%%reactor. Please ensure that there are as many transformers and reactors
%%as specified in the main. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function initialise_systemdata(topology)
global CONSTANTS mpopt Systemdata Optimisation

%%Surpress MATPOWER outputs
mpopt = mpoption('verbose',0,'out.all',0);

%%Structure containing power system and optimization related information
Systemdata.mpc = topology;
Systemdata.Nbranch = size(Systemdata.mpc.branch,1);
Systemdata.Nbus = size(Systemdata.mpc.bus,1);
Systemdata.Nstring = size(Systemdata.mpc.gen,1);

%%Logic vector with 1 on wtg positions in gen matrix
Systemdata.wtg_pos = logical([0; ones(Optimisation.Nturbines,1); ...
    zeros(Optimisation.Npv,1)]); 

%%If Nturbines = 91, the string model is used. Then, the turbine types are
%%extracted from the excel file
if Optimisation.Nturbines == 91
    [~, types] = xlsread('Branch_calculations_turbinelevel.xlsx','Sheet1', 'O146:O251');
    Systemdata.wtg_type = types(not(strcmp('',types)));
end                                          
%%logic vector with 1 on pvg positions in gen matrix
Systemdata.pvg_pos = logical([zeros(Optimisation.Nturbines+1,1); ...
    ones(Optimisation.Npv,1)]);

%%Compute the transformer lookup table if the transformers are controlled
if Optimisation.Ntr ~=0
    %%Logic vector with 1 on transformer positions in matrix 
    %%Finc transformers
    Systemdata.trans = Systemdata.mpc.branch(:,CONSTANTS.ANGMAX) ~= 0;
    
    %%Extract the minimum and maximum ratio's
    R(:,1) = Systemdata.mpc.branch(Systemdata.trans,CONSTANTS.ANGMAX);
    R(:,2) = Systemdata.mpc.branch(Systemdata.trans,CONSTANTS.ANGMIN);
    
    %%Compute a lookup table for the transformers for the number of taps
    %%per transformer
    Systemdata.trlookup = NaN * ones(Optimisation.Ntr,max(Optimisation.Ntaps));
    for i = 1:Optimisation.Ntr
        Systemdata.trlookup(i,1:Optimisation.Ntaps(i,1)) = linspace(R(i,1),R(i,2),Optimisation.Ntaps(i,1));
        Systemdata.trstep(i) = abs(Systemdata.trlookup(i,1)-Systemdata.trlookup(i,2));
    end
else
    Systemdata.trans = [];
    Systemdata.trstep = [];
end

%%Compute the positions of the shunt reactors only if the shunts are
%%controlled
if Optimisation.Nr ~=0
    %%Logic vectors with 1 on reactor positions in bus and branch matrices
    %%Find the busses of the reactors
    Systemdata.shunts = Systemdata.mpc.bus(:,CONSTANTS.BS)~= 0;
    
    %%Find the branches to which the reactors belong
    Systemdata.shuntbranch = (Systemdata.mpc.branch(:,CONSTANTS.T_BUS)...
    == Systemdata.mpc.bus(Systemdata.shunts,CONSTANTS.BUS_I)|Systemdata.mpc.branch(:,CONSTANTS.F_BUS)...
    == Systemdata.mpc.bus(Systemdata.shunts,CONSTANTS.BUS_I))~=0;
else
    Systemdata.shunts = [];
    Systemdata.shuntbranch=[];
end
end
