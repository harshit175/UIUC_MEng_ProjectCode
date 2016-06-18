global n1 n2 n3 n m totalnodes

i=1:totalnodes;
j=totalnodes+1:2*totalnodes;

add_to_te=10000;

T(totalnodes,te/recordstep+1:add_to_te/recordstep)=0;
alpha(totalnodes,te/recordstep+1:add_to_te/recordstep)=0;



starttime=cputime;
for step=te+1:te+add_to_te
  tic    
Pc_pre=Pcfunction(U(i),U(j),mlump,1);
    
    if beta~=0
    disp('entering implicit algorithm')
    implicitscheme %uncoupled
    disp('exiting implicit algorithm')
    U(i)=Uit(i,iteration); U(j)=Uit(j,iteration);
    else        
    U = A\(B*U)+A\sparse(vertcat(F', Pc_pre));
    end
    if not(isreal(U))
        disp('curing values imaginary at time step');
        disp(step)
        complexcurestep=step;
        break
    end
    %T(:,step)=U(i);
    %alpha(:,step)=U(j);
        
    %if check==0 && (exist('bb'))
    %Pc(bb)=0;
    %end
    %the implicit scheme will be coded here, it will have U as input
    %and Pc_implicit as output
    
    
    if mod(step,recordstep)==0
        T(:,step/recordstep)=U(i);
        alpha(:,step/recordstep)=U(j);
        time(step/recordstep)=step*dt;
        if (mean(alpha(nodesnearignition,step/recordstep))>curelimit && check==0)
            A=A0; A11=A011; B=B0; B11=B011; Im1=Im01;  %these conversions not needed for ignition by heating
            F=F0; 
            check=1;
            disp('ignition stopped at time step');
            disp(step)
            ignoffstp=step;
        end
        
        if(mean(alpha(nodesnearignition,step/recordstep))>curelimit && mean(alpha(nodesnearresist,step/recordstep))>curelimit)
        F=F0;
        elseif(mean(alpha(nodesnearignition,step/recordstep))>curelimit && mean(alpha(nodesnearresist,step/recordstep))<curelimit)
            F=Fresistive;
        elseif(mean(alpha(nodesnearignition,step/recordstep))<curelimit && mean(alpha(nodesnearresist,step/recordstep))>curelimit)
            F=Fignition;
        end
    end
    
    if mod(step,recordstep/10)==0        
        step
    end
   toc 
end
endtime=cputime-starttime;
te=step;