%Post processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tnode=zeros(m,te/recordstep);
dTnode=zeros(m,te/recordstep);
for i=1:te/recordstep
for sensor=1:1+n1;
    Tnode(sensor,i)=Tm(3,sensor+n3,i);
   
end
end

for i=2:te/recordstep
for sensor=1:1+n1;
    dTnode(sensor,i)=(Tnode(sensor,i)-Tnode(sensor,i-1))/(dt*recordstep);   
end
end

[dTnodemax,Fronttimestep]=max(dTnode');

Axisdistanceset=1e3*(rin:deltar(1):(rinter));
Fronttime=horzcat(dt*recordstep*Fronttimestep);

plot(Fronttime,Axisdistanceset, 'k*')

xlabel('Time (sec.)')
ylabel('Front position (mm)')
title ('Front position vs time for r= mm')

% s1=15;  s2=35;
% for i=1:te
% %Tnode1(i)=Tm(1,1,i);
% Tnode1(i)=Tm(s1,1,i);
% Tnode2(i)=Tm(s2,1,i);
% alphanode1(i)=alpham(s1,1,i);
% end
% 
% Tdiffnode1=0;
% Tdiffnode2=0;
% for i=1:te-1
% Tdiffnode1(i+1)=(Tnode1(i+1)-Tnode1(i))/dt;
% Tdiffnode2(i+1)=(Tnode2(i+1)-Tnode2(i))/dt;
% end

% for i=1:te
% time(i)=dt*(i-1);
% end

%Ttimepoint=Tm(1:m,2,500);

%surface plot of unsteady state.
%unsteadyplotanimated
%aaa
