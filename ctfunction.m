function [ ct ] = ctfunction(rleft,rright,zbot,ztop,mat,dt)
%This function calculates the CC as given in the Qi Zhu paper and then mass
%lumps it.

j=0.25*(rright-rleft)*(ztop-zbot);
%j=1;
ct=j^2*1/dt*[ (4*(3*rleft + rright))/(9*(rleft - rright)*(zbot - ztop)),   (4*(rleft + rright))/(9*(rleft - rright)*(zbot - ztop)),   (2*(rleft + rright))/(9*(rleft - rright)*(zbot - ztop)), (2*(3*rleft + rright))/(9*(rleft - rright)*(zbot - ztop));
   (4*(rleft + rright))/(9*(rleft - rright)*(zbot - ztop)), (4*(rleft + 3*rright))/(9*(rleft - rright)*(zbot - ztop)), (2*(rleft + 3*rright))/(9*(rleft - rright)*(zbot - ztop)),   (2*(rleft + rright))/(9*(rleft - rright)*(zbot - ztop));
   (2*(rleft + rright))/(9*(rleft - rright)*(zbot - ztop)), (2*(rleft + 3*rright))/(9*(rleft - rright)*(zbot - ztop)), (4*(rleft + 3*rright))/(9*(rleft - rright)*(zbot - ztop)),   (4*(rleft + rright))/(9*(rleft - rright)*(zbot - ztop));
 (2*(3*rleft + rright))/(9*(rleft - rright)*(zbot - ztop)),   (2*(rleft + rright))/(9*(rleft - rright)*(zbot - ztop)),   (4*(rleft + rright))/(9*(rleft - rright)*(zbot - ztop)), (4*(3*rleft + rright))/(9*(rleft - rright)*(zbot - ztop))];


%if mat==1
   ct=sum(sum(ct))/sum(diag(ct))*[ct(1,1) 0 0 0; 0 ct(2,2) 0 0; 0 0 ct(3,3) 0; 0 0 0 ct(4,4)]; 
%end
end

