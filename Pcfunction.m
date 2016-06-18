
function [ Pc ] = Pcfunction(Tstep,alphastep,mlump, para )
%PCFUNCTION Summary of this function goes here
%   Detailed explanation goes here

%this function implements the cure kinetics formula, and is the only
%function that will need to be changed if a different chemistry is used

 %para means parameter
 
global totalnodes n1 n3 n m bb multfactor


Pc=zeros(totalnodes,1);
ec=exp(1);
pc=0;

%k01=ec^20.2;
%k02=ec^9.5;

k01=600*ec^20.2; %kinetics accelerated by 600 times to match the experimental front speed
k02=600*ec^9.5;

Eatt1=64735;
Eatt2=31561;

mm=1;
nn=1.07;

p=6.96e-3;
q=4.57e-1;

R=8.314;

lumpingcounter=0; 
for c=n3+1:n1+n3+1 
    %only radial nodes in the channel are considered
     
%     if (c>=n3+1 && c<=n1+n3)  %will have to be modified
%         mat=1; channel=1;
%     else
%         mat=2; channel=0;
%     end
    for r=1:m+1 % moving along axial direction
        
        
        a=(r-1)*(n+1)+c;  %this is the bottom left node number related to element number
        
        if r==1 %increment done only when r==1, tis way the structuredness of the mesh gets utilized
            
            lumpingcounter=lumpingcounter+1;
            %mlumpelement=mlump(:,lumpingcounter);
            %     r
            %     c
            
        end
        

        
        
        for i=1:1
            k1=k01*ec^(-Eatt1/(R*(Tstep(a)+273)));
            k2=k02*ec^(-Eatt2/(R*(Tstep(a)+273)));
            
            dk1bydT=k1*(Eatt1/R)*(1/(Tstep(a)+273))^2;
            dk2bydT=k2*(Eatt2/R)*(1/(Tstep(a)+273))^2;
            
            if Tstep(a)<78
                alphmax=p*Tstep(a)+q;
                if alphmax-alphastep(a)<0
                    alphmax=alphastep(a);
                    %disp('curing value had to be adjusted')
                    %alphastep(a)
                    %Tstep(a)
                    %disp('row number')
                    %r
                    %disp('column number')
                    %c
                end
            else
                alphmax=1;
            end
            if para==1
            %pc(i)=mlump(i)*(1/60)*(k1+k2*(alphastep(a))^mm)*(alphmax-alphastep(a))^nn;
            pc(i)=(1/60)*(k1+k2*(alphastep(a))^mm)*(alphmax-alphastep(a))^nn; % PC expression
            
            elseif para==2
                pc(i)=(1/60)*(dk1bydT+dk2bydT*(alphastep(a))^mm)*(alphmax-alphastep(a))^nn; %derivative of Pc w.r.t. temperature
                %pc(i)=dpcbydT(i);
            
            
            else
                pc(i)=(1/60)*((k2*mm*alphastep(a)^(mm-1))*(alphmax-alphastep(a))^nn - nn*(k1+k2*(alphastep(a))^mm)*(alphmax-alphastep(a))^(nn-1));
                %pc(i)=dpcbydalpha(i); %derivative of Pc w.r.t alpha
            end
            % if not(isreal(pc(i)))
            %     Disp('Improve definition of Pcfunction')
            % end
        end
        
        if (r==1 || r==m+1)
            Pc(a)=multfactor(lumpingcounter)*pc; 
        else
            Pc(a)=2*multfactor(lumpingcounter)*pc; %multiplying factor is twice in inner row, compared to top and bottom row
        end
        
        %Pc(bb)=0; %This line when activated, fixes the cure at bottom nodes to the initial value (usually zero)     
    end
end





end
%lumpingcounter


