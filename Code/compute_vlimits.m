%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This funciton calculates the voltage limits at the slackbus as function
%%of the reactive power at the PCC. For this, the requirement of the grid 
%%code is used. If the reactive power is outside the specified range, the
%%voltage gets the value specified at the bounds (Qpcc = +/- 0.99)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vlimits = compute_vlimits(Qpcc) 
%%Check bounds of Qpcc
if Qpcc > 1
    Qpcc = 0.99;
elseif Qpcc < -1
    Qpcc = -0.99;
end

%%Corner coordinates of QV region
x = [-0.25 1 1 0 -1 -1];
y = [0.9 0.9 1 1.1 1.1 1];

%%Make shape
shape = polyshape(x,y);

%%Add a vertical line specified by Qpcc
line = [Qpcc 0.8; Qpcc 1.2];

%%Compute the intersection of the region and the vertical line i.e. the
%%allowed bounds of v at the specified Q
[int] = intersect(shape,line);
vlimits = flip(transpose(int(:,2)));
end