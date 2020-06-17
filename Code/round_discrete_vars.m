%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to round the discrete variables of 'Xin:
%%-Transformer positions in the solution vector are rounded using
%%lookup tables describing all possible ratios
%%-Reactors are rounded to 1/0 depending on on/off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xout = round_discrete_vars(Xin)
global Optimisation Systemdata;

%%Copy the output to the input. This is for the continuous variables.
Xout = Xin;

%%Round the discrete variables. For this, the functions checks per variable
%%of Xin whether it is discrete or not. 
for i = 1:Optimisation.Nvars
    %%If the variable is a disctrete transformer, round it
    if Optimisation.discrete(i) & Optimisation.tr_pos(i)
        %%Offset used to calculate the entry in the transformer tap
        %%positions lookup table
        tri = i - Optimisation.Nturbines - Optimisation.Npv;
        
        %%Round the transformer tap ratios to the nearest tap i.e. the
        %%nearest ratio from the transformer lookup table.
        Xout(i) = interp1(Systemdata.trlookup(tri,1:Optimisation.Ntaps(tri,1)),...
        Systemdata.trlookup(tri,1:Optimisation.Ntaps(tri,1)),Xin(i),'nearest');
    
    %%If the variable is a discrete reactor, round it
    elseif Optimisation.discrete(i)& Optimisation.r_pos(i)
        Xout(i) = round(Xin(i));
    end
end    
end