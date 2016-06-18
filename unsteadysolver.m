%Data from input file will be used. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Processing starts

%The idea here is that, PDMS has element length increasing with geometric progression in radial direction, with
%ratio equal to gpratio
if gpratio~=1
    n2int=int8(log(1-(rout-rinter)/(rinter-rin)*(1-gpratio)*n1)/log(gpratio));
    %calculating number of elements
    n2int=n2int+1;
    %increasing elements by one
    n2=double(n2int);
    
    for i=1:n2+1
        r2(i)=rinter+ (rinter-rin)/n1*(1-gpratio^(i-1))/(1-gpratio); %i=1, r2=rinter; i=2, r2=rinter+a;  i=3, r2=rinter+a+a*gpratio and so on
        deltar=[(rinter-rin)/n1 (rout-rinter)/n2 (rin-rc)/n3];
    end
    
else
    n2=double(int64((rout-rinter)*1e3*000));  %n2 as to be specified here if gpratio =1, the code will be erroneous otherwise,
    %if n2 ==0 with gpratio 1, then PDMS won't be a part of the domain
    if n2==0 %n2 can be kept at zero if gpratio is not equal to 1, 
        deltar=[(rinter-rin)/n1 0 (rin-rc)/n3];
    else
        deltar=[(rinter-rin)/n1 (rout-rinter)/n2 (rin-rc)/n3];
    end
end
%if PDMS is not to be a part of the domain, keep gpratio as 1 and n2 as
%zero

n=n1+n2+n3; %total number of radial elements


global totalnodes;
totalnodes =(m+1)*(n+1);
curingnodes=(m+1)*(n1+1);



rinm=[rin rinter rc]; %in stands for inner radius, m stands for material

%For instance, innr radius for PDMS (mat=2) is rinter
%the order is [microchannel PDMS conductor] (not the order in which they are placed radially)
%the code was initially developed without the wire, hence the discrepancy
%in ordering

deltaz=hz/m; % element length in axial direction

%Assembling KT and CT matrix
K=sparse(totalnodes,totalnodes);
CT=sparse(totalnodes,totalnodes);
Cc=sparse(totalnodes,totalnodes);



F(totalnodes)=0; %load vector initialised at zero, this will be modified later and used in the
%ignition stage
F0=F;  %this will be used after ignition boundary is insulated


disp('matrix assembly starting')
i=1:4; lumpingcounter=1;
mlump=[0; 0; 0; 0];
for c=1:n
    if (c<=n3)
        mat=3; channel=0; %conductor
    elseif (c>=n3+1 && c<=n3+n1)
        mat=1; channel=1; %channel
    else
        mat=2; channel=0; %PDMS
    end
for r=1:m
    if mat==1
        rleft=rinm(mat)+(c-n3-1)*deltar(mat);
        rright=rleft+deltar(mat);
    elseif mat==2 
        if gpratio~=1 %i.e. PDMS, has two cases, one for a gp in the r direction, and other is no gp
            rleft=r2(c-n1-n3);
            rright=r2(c-n1-n3+1);
        else
            rleft=rinm(mat)+(c-n1-n3-1)*deltar(mat);
            rright=rleft+deltar(mat);
        end
    elseif mat==3 %conductor
        rleft=rinm(mat)+(c-1)*deltar(mat);
        rright=rleft+deltar(mat);
    end
    
    zbot=(r-1)*deltaz;
    ztop=r*deltaz;
    
    if r==1 %this if condition to be used only when mesh in z direction is uniform
    k=stiffnessfunction(rleft,rright,zbot,ztop,mat,kc); %stiffness matrix for a single element
    ct=rhocp(mat)*ctfunction(rleft,rright,zbot,ztop,mat,dt);    %ct for a single element
    end
    
    a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c+1 r*(n+1)+c]; % nodes: 1, 2, 3, 4 (traversed anticlockwise)
    
    
    %r*(n+1)+c   --------------   r*(n+1)+c+1
    %    |                              |
    %    |                              |
    %    |         Element(r,c)         |
    %    |                              |
    %    |                              |
    %(r-1)*(n+1)+c -------------  (r-1)*(n+1)+c+1
   
    %local to global assembly of k and ct
    K(a(i),a(i))=K(a(i),a(i))+k(i,i); 
    CT(a(i),a(i))=CT(a(i),a(i))+ct(i,i);
    
   %local to global assembly of cc
    cc=channel*ctfunction(rleft,rright,zbot,ztop,mat,dt);
    Cc(a(i),a(i))=Cc(a(i),a(i))+cc(i,i);
    
    if mat==1 && r==1  %capitalizes on the structured nature of microchannel mesh, so is only computed for the first row of elements.
        mlump(:,lumpingcounter)=diag(ctfunction(rleft,rright,zbot,ztop,mat,1)); %mlump is a 4 x (lumpingcounter) matrix
        lumpingcounter=lumpingcounter+1;  %lumping counter's final value is equal to the radial nodes in the microchannel       
    end
    
    if mat==3 %this function is to add resistive heating to the conductive wire
        qhe=heatsourcevectorfunction( rleft,rright,zbot,ztop,Q);
        F(a(i))=F(a(i))+qhe(i)';
    end
    
end
c/n %this displays the progress of assembly as a fraction
end

Fresistive=F;

 %calculation of multiplication factors (multfactor) for Pcfunction
%multfactor is basically the factor disussed in eq. 26 of the Zhu paper
%the factors are calculated for a single row of nodes


%the algorithm is designed to speed up the computation of Pc in Pcfunction
%and exploits the structured nature of the mesh to a great deal

lumpingcounter=0; 
for cc=n3+1:n1+n3+1
    lumpingcounter=lumpingcounter+1;
    if (cc==n3+1 ) %leftmost radial node in channel
        multfactor(cc-n3)=mlump(1,lumpingcounter); %leftmost node
        
    elseif(cc==n3+n1+1) %rightmost radial node in channel
        multfactor(cc-n3)=mlump(2,lumpingcounter-1); %rightmost node
        
    else %innder radial nodes in channel, each node is common between two elements, so gets contribution from left and right
        multfactor(cc-n3)=mlump(1,lumpingcounter)+mlump(2,lumpingcounter-1); %innder nodes, get contribution from left and right
        
    end
end

%The zero diagonal elements have to be changed to 1, else the Cc matrix is
%singular, this does not change the solution as the alpha corresonding to
%those elements is in PDMS (mat=2)or copper (mat=3) and is a constant (equal to zero)
for i=1:totalnodes
    if (Cc(i,i)==0)
        Cc(i,i)=1;
    end
end
%Assembly of KT, CT and Cc done
disp('assembly of matrices complete')
%initialising alpha to a small value<<<1 at t=0
% i=1:4;
% for c=1:n
%     if (c>=n3+1 && c<=(n1+n3))
%         mat=1; curingnode=1;
%     elseif (c>=n1+n3+1)
%         mat=2; curingnode=0;
%     elseif (c<=n3)
%         mat=3; curingnode=0;
%     end
%     
% for r=1:m
%     a=[(r-1)*(n+1)+c (r-1)*(n+1)+c+1 r*(n+1)+c r*(n+1)+c+1];
%     alpha(a(i),1)=curingnode*0.000;
% end
% end
%alpha=alpha';
%initialization complete






%Assembling load vector (to ignite the monomer), 
%simulating soldering iron 
i=1:4;
for c=1:(n3+n1)
    if (c<=n3)
        mat=3; %conductor
        rleft=rinm(mat)+(c-1)*deltar(mat); %calculates rleft
    elseif (c>=n3+1 && c<=n3+n1)
        mat=1; channel=1; %channel
        rleft=rinm(mat)+(c-n3-1)*deltar(mat); 
    end
    
    rright=rleft+deltar(mat); %calculates rright, deltar(mat) has two values, one for each material
    Qve=heatvectorfunction(rleft,rright,q);
    
    a=[c c+1 c+2 c+3];  
    F(a(i))=F(a(i))+Qve(i);
    
    
end

Fresistiveplusignition=F; %this is the load vector when there is both resistive heatig and ignition
Fignition=Fresistiveplusignition-Fresistive; %load vectors for ignition only
%F(c+1)=0;
%Load vector assembled


%Applying constant temperature boundary conditions on right boundary

j=1:(m+1); %Axial node set

global bb
bb=1:n1+n3+1;   %bb=ignition boundary nodes


    
    


br=(n+1)*j; %br=right boundary nodes
bl=br-n;    %bl=left boundary nodes




nodesnearresist=horzcat(bl+n3,bl+n3+1,bl+n3+2); %the column of nodes near the wire
nodesnearignition=n+2+n3:n+2+n1+n3; %the row of nodes one above the ignition boundary


Pc_pre=zeros(totalnodes,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solving for unsteady state

%implicit solver

%The idea is to have two sets of matrices, one for post ignition, and one
%for ignition

%Matrices for post ignition have all DOFs free, so they have a subscript 0
%0 (nought)

A11=sparse(CT+beta*K);  %has to be modified for ignition

A12=sparse(-rho1*Hr*Cc);

A21=sparse(totalnodes,totalnodes);

A22=sparse(Cc);

A=sparse([A11 A12; A21 A22]);  %has to be modified for ignition


A0=A;        %stored, since it will be modified for ignition
A011=A11;    %stored, since it will be modified for ignition


B11=sparse(CT-(1-beta)*K);  %has to be modified for pre ignition

B12=sparse(-rho1*Hr*Cc);

B21=sparse(totalnodes,totalnodes);

B22=sparse(Cc);

B=[B11 B12; B21 B22];  %has to be modified for ignition

B011=B11;%stored, since B11 will be modified for ignition
B0=B;%stored, since B will be modified for ignition






T(totalnodes,te/recordstep)=0; alpha(totalnodes,te/recordstep)=0; %T and alpha vectors (matrices if time steps considered) initialized.


T(1:totalnodes,1) = Tambient; %inital temperature condition

alpha(1:totalnodes,1) = alphainitial; %initial curing condition

Tinitialrows=double(int64(m/10)); %This is a relic,  to prove that ignition temperature BC was causing oscillations

for count=1:Tinitialrows %this is for a graded initial condition
   %T(count*(n+1)+1:count*(n+1)+n3+1)=Tambient+(Tinitialrows-count)*(Tb-Tambient)/Tinitialrows; 
end

%bb stands for bottom boundary
%note that the code deletes bb (variable containing bottom bounary node numbers)
%if there is ignition (i.e. q is not zero)
if q~=0  
   clear bb 
end 

if (exist('bb'))  
    for p=1:totalnodes
        F(p)=F(p)+Tb*( sum(B11(p,bb)) - sum(A11(p,bb)) ); %needs to be commented out if ignition 
        %is to be done using heating and not temperature bc
        
    end
    
    F(bb)=0; %the forcing at constant temperature boundary is made zero
    
    T(bb,1)=Tb;  %sets the temperature of bottom boundary
    
    
    %nodesnearignition=n+2+n3;
        
    A11=modifymatrices(A11,bb); %this gives A11 for ignition BC
    %A12=modifymatrices(A12,bb);
    %A22=modifymatrices(A22,bb);
    
    A=sparse([A11 A12; A21 A22]); %this gives A for ignition BC
    
    %invA=inv(A);
    %invA = A\IA;
    B11=modifymatrices(B11,bb);
    %B12=modifymatrices(B12,bb);
    %B22=modifymatrices(B22,bb);
    
    B=[B11 B12; B21 B22]; %this gives B for ignition BC
    
    %invAB=invA*B;
end

i=1:totalnodes;
j=totalnodes+1:2*totalnodes;  
%these two indices separate the U vector into temperature and cure vectors
%ignitionnodes=1:n1;



Im1=sparse([A11 A12]); %during ignition, gets passed to implicitscheme file
Im01=sparse([A011 A12]); %post ignition, gets passed to implicitscheme file

U=vertcat(T(:,1), alpha(:,1)); %U is the composite vector of temperature and cure vectors

check=0; %this gets changed to 1 after ignition is switched off
%Pc=zeros(totalnodes,1);
disp('time marching starting')

for step=2:te
tic    
Pc_pre=Pcfunction(U(i),U(j),mlump,1); %calculates the forcing for curing, i.e. the cure heat source term
    
    if beta~=0
    disp('entering implicit algorithm')
    implicitscheme %this applies the Newton Raphson iteration
    disp('exiting implicit algorithm')
    U(i)=Uit(i,iteration); U(j)=Uit(j,iteration);
    else        
    U = A\(B*U)+A\sparse(vertcat(F', Pc_pre)); %the time consuming inversion step in the explicit scheme
    end
    if not(isreal(U)) %complex cure values break the time stepping
        disp('curing values imaginary at time step');
        disp(step)
        Ucomplex=U;
        complexcurestep=step; %the step at which values go complex gets recorded
        break
    end
    %T(:,step)=U(i);
    %alpha(:,step)=U(j);
        
       
    if mod(step,recordstep)==0
        %U is  broekn down into T and alpha 
        T(:,step/recordstep)=U(i); %Values are recorded only at 'recordstep' intervals 
        alpha(:,step/recordstep)=U(j);
        time(step/recordstep)=step*dt; %the time vector also gets populated as time marching progresses
        %the following loop checks cure value near ignition and switches
        %off the ignition
        if (mean(alpha(nodesnearignition,step/recordstep))>curelimit && check==0)
            A=A0; A11=A011; B=B0; B11=B011; Im1=Im01;  %these conversions not needed for ignition by heating
            F=F0; 
            check=1; %this stops subsequent execution of the if condition
            disp('ignition stopped at time step');
            disp(step)
            ignoffstp=step; %the step at which ignition stops gets recorded
        end
        
        %the below 3 conditions ensure that the right forcing function is
        %applied
        if(mean(alpha(nodesnearignition,step/recordstep))>curelimit && mean(alpha(nodesnearresist,step/recordstep))>curelimit)
        F=F0;
        
        elseif(mean(alpha(nodesnearignition,step/recordstep))>curelimit && mean(alpha(nodesnearresist,step/recordstep))<curelimit)
            F=Fresistive;
            
        elseif(mean(alpha(nodesnearignition,step/recordstep))<curelimit && mean(alpha(nodesnearresist,step/recordstep))>curelimit)
            F=Fignition;
            
        end
    end
    
    if mod(step,recordstep/10)==0        
        step %This prints out the step at certain intervala 
    end
   toc 
end



%Deleting variables not needed any more

%clear check i j
%clear K CT Cc DT 
%clear K0 CT0 Cc0 DT0 
%clear A B AB A0 B0 AB0

%clear A A0 A011 A11 A12 A21 A22 B B0 B11 B12 B21 B22 Im1 Im01 Im1 Im21 Im22
%clearing variables no longer needed saves memory

Axisdistanceset=1e3*(0:deltaz:hz);

%postprocessing in 'postprocessor' file







