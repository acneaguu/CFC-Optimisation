%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to compute the allowed deviation of Qpcc w.r.t. the
%%Qsetpoint requested by the TSO. In this function, it is assumed that the
%%gridcode of Qpcc is rectangular. Otherwise, Qmin and Qmax need to be
%%computed by taking the intersect of a shape (corresponding to the grid
%%code) and the active power output. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function qpcc_limits(Qsetpoint)
global Qref

%%Grid code for rectangle
Qmin = -0.4;
Qmax = 0.33;

%%Compute the limits
Qref.limits = [max([Qmin, Qsetpoint-Qref.tolerance]), ...
    min([Qmax,Qsetpoint+Qref.tolerance])];
end