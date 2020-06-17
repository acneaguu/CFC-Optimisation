%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to compute the active and reactive power output of
%%Nstrings for a given solar irradiance. The output is a vector containing
%%the same power output for the different strings. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ptot,Qtot] = compute_pq_pvg(irradiance,Nstrings)
%%Parameters
%Module specs
area_module = 1.046*1.558;      %m2
effiency_module = 0.221;
Nmodules_per_string = 41667;    
Prated = 360;                   %W

%Inverter specs
efficiency_inverter = 0.95;
Qrange = 1/3;                   %percentage of Pavailable

%%Compute P for a single module
P = irradiance * effiency_module * area_module;               
if P > Prated
    P = Prated;
end

%%Compute P and Q total for a string
Ptot = P * Nmodules_per_string* efficiency_inverter * 1e-6; %P in MW 
Qtot = Qrange * Ptot;                                       %Q in MVAr

%%Repeat for N strings
Ptot = repmat(Ptot,1,Nstrings);
Qtot = repmat(Qtot,1,Nstrings);
end