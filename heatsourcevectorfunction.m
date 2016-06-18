function [ qhe ] = heatsourcevectorfunction( rleft,rright,zbot,ztop,Q)
%IF the heating is done resistively in the wire, then function is used by the code, it
%takes the heating Q in W/m^3 as input and outputs the the heat source
%vector for the 4 nodes in an element.

qhe=Q*[ ((rleft - rright)*(2*rleft + rright)*(zbot - ztop))/12, 
    ((rleft - rright)*(rleft + 2*rright)*(zbot - ztop))/12, 
    ((rleft - rright)*(rleft + 2*rright)*(zbot - ztop))/12, 
    ((rleft - rright)*(2*rleft + rright)*(zbot - ztop))/12];


end

