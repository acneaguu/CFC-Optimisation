%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to save the workspace. Data is an optional input
%%variable which is saved. It is possible to specify a prefix which is
%%saved before the rest of the file name. The format is as follows:
%%userstring_With/No Opt_Algorithm name_Number of runs_Number of
%%variables_date + time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function savedata(Data)
global CONSTANTS Qref mpopt Systemdata PFresults Optimisation Results Keeptrack FCount;

%%Name of the struct  
%name = input('Enter custom prefix:\n','s');
name = 'system_13_paramsweep';

%%Check whether there was optimisation or not
if or(or(Optimisation.w1,Optimisation.w2),Optimisation.w3)
    optstr = 'With Opt ';
else
    optstr = 'No Opt ';
end

%%Check which algorithm
if Optimisation.algorithm == 1
    algstr = 'GA_';
elseif Optimisation.algorithm == 2
    algstr = 'PS_';
elseif Optimisation.algorithm == 3
    algstr = 'CDEEPSO';
elseif Optimisation.algorithm == 4
    global parameter
    configmvmo = sprintf('fs=%4.1d and %4.1d_Nmut=%4.1d and %4.1d',...
        parameter.fs_factor_start,parameter.fs_factor_end,...
        parameter.n_random_ini,parameter.n_random_last);
    algstr = strcat('MVMO-SHM',configmvmo);
end

%%Create an empty string if name is nonexistent
if isempty(name)
    name = [];
else
    name = strcat(name,"_");
end

%%Make a string containing some variables of a run
rundata = sprintf('Nruns =%3.1d Nvars =%3.1d Nswarm =%3.1d',...
    Optimisation.Nruns,Optimisation.Nvars,Optimisation.Populationsize);

%%Makes the final name of the data and saves it
namestr = strcat(name,optstr," ",algstr,rundata,"_",datestr(now,'dd-MM-yyyy HH-mm-ss'),'loc');
save(namestr)
end