%Post processing, contains front velocity code as well
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tnode=zeros(m,te/recordstep);
dTnode=zeros(m,te/recordstep);
for i=1:te/recordstep
for sensor=1:m+1;
    radialnode = 26; %the node from the centre
    Tnode(sensor,i)=Tm(sensor,radialnode,i);
   
end
end

for i=2:te/recordstep
for sensor=1:m+1;
    dTnode(sensor,i)=(Tnode(sensor,i)-Tnode(sensor,i-1))/(dt*recordstep);   
    %this is dT/dt at every node along an axial line passing through radialnode
end
end

[dTnodemax,Fronttimestep]=max(dTnode'); 
%dTnodemax gives the distance of Front
%Fronttimestep gives the time when the front is at a certain distance

Axisdistanceset=1e3*(0:deltaz:hz);
Fronttime=horzcat(dt*recordstep*Fronttimestep);
for(frontcount=1:length(Fronttime)) %tis loop ensures that front distance vs time plot 
   %is only plotted till front distance and not entire length
    if (abs(Fronttime(frontcount)-dt*step)<1e-10)
        frontcount
        break;
    end
end
back=6;
for count=back+1:frontcount
FrontVel(count)=(Axisdistanceset(count)-Axisdistanceset(count-back))/(Fronttime(count)-Fronttime(count-back));
end

%smoothing front velocity
smoothFrontVel=smooth(smooth(smooth(FrontVel)));
figure(1)
plot(Fronttime(1:count),Axisdistanceset(1:count))
xlabel('Time (sec.)','fontsize',14)
ylabel('Front position (mm)','fontsize',14)
%title ('Front position vs time for r= mm')

figure(2)
plot(Fronttime(1:count),smoothFrontVel(1:count))
xlabel('Time (sec.)','fontsize',14)
ylabel('Front velocity (mm/s)','fontsize',14)

% ft=Fronttime(1:count);
% p1=-9.4;
% p2=30;
% p3=-37;
% p4=22;
% p5=-4.3;
% p6=1.4;
% FrontVel=p1*6*ft.^5+p2*5*ft.^4+p3*4*ft.^3+p4*3*ft.^2+p5*2*ft.^1+p6;


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