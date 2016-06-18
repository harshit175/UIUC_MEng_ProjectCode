function [ Qve ] = heatvectorfunction( rleft,rright,Q )
%This function adds the forcing function for ignition at bottom boundary,
%since only two nodes out of the four elements in a node are ignited, the
%last two entries are zero.


Qve=Q*[ -((rleft - rright)*(2*rleft + rright))/6, -((rleft - rright)*(rleft + 2*rright))/6, 0, 0];
end

