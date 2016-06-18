%the code implements the newton raphson scheme
clear Uit

Uit(:,1)=(U); %Uit means U that is iterated


P_pre=sparse(vertcat(F', Pc_pre));
gamma=B*U+(1-beta)*P_pre;
clear Residue

for iteration=1:1000
    iteration
    
    Pc(:,iteration)=Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,1);
    Pit=(vertcat(F', Pc(:,iteration))); %it means iteration
    
    
    Residue(:,iteration)=A*Uit(:,iteration)-beta*Pit-gamma; 
    
    if iteration==1
        firstitnorm=norm(Residue(i,1));
    end
    if norm(Residue(i,iteration))<.0001*firstitnorm
        break
    end
    %Im21=-sparse(beta*diag(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,2)));
    Im21=-beta*sparse(1:totalnodes,1:totalnodes,(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,2))); 
    %Im21 is the derivative of d(alpha)/dt w.r.t. T (off diagonal
    %jacobian)
    if check==0 && exist('bb')
        %Im21=modifymatrices(Im21,bb); %This step fixes the curing at ignition nodes
    end
    
    %Im22=sparse(Cc-beta*diag(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,3)));
    Im22=Cc-beta*sparse(1:totalnodes,1:totalnodes,(Pcfunction(Uit(i,iteration),Uit(j,iteration),mlump,3)));
    %Im22 is the derivative of d(alpha)/dt w.r.t alpha (Diagonal jacobian)
    
    %tt2=cputime;
    Im=[Im1; Im21 Im22];
    %timeelapsed2=cputime-tt2
    
    %tic
    %disp('dU calculation starts')
    dU=-(Im\Residue(:,iteration));  %this inversion step is the most time consuming step in the whole process !!
    %disp('dU calculation ends')
    %toc
    
    
    %dUrecord(:,iteration,step)=dU;
    Uit(:,iteration+1)=Uit(:,iteration)+dU;    
end %this is the Newton-Raphson iteration loop






