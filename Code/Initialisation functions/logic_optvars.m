%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function creates multiple logical vectors which indicate the
%%position of the turbines, pv strings, transformer positions and reactor
%%positions. Moreover, it makes a vector indicating which variables are
%%discontinuous. This is used to round the discrete variables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function logic_optvars()
global Optimisation

%%Logic vector with 1 on wtg positions
Optimisation.wtg_pos = logical([ones(1,Optimisation.Nturbines) ...
    zeros(1,(Optimisation.Npv+Optimisation.Ntr+Optimisation.Nr))]); 

%%Logic vector with 1 on transformer positions                                             
Optimisation.pvg_pos = logical([zeros(1,Optimisation.Nturbines) ...
    ones(1,Optimisation.Npv) zeros(1,(Optimisation.Ntr+Optimisation.Nr))]); 

%%Logic vector with 1 on transformer positions                                         
Optimisation.tr_pos = logical([zeros(1,(Optimisation.Nturbines +...
    Optimisation.Npv)) ones(1,(Optimisation.Ntr)) zeros(1,Optimisation.Nr)]); 

%%Logic vector with 1 on reactor positions
Optimisation.r_pos = logical([zeros(1,(Optimisation.Nturbines + ...
    Optimisation.Npv + Optimisation.Ntr)) ones(1,Optimisation.Nr)]); 

%%Logic vector which is 1 for discrete variables
Optimisation.discrete = zeros(1,Optimisation.Nvars);
Optimisation.discrete(Optimisation.which_discrete) = 1;
Optimisation.discrete = logical(Optimisation.discrete);
end