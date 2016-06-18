clear count
for count=1:te/recordstep
    Tm(:,:,count)=zeros(m+1,n+1);  %Temperature matrix initialized
    alpham(:,:,count)=zeros(m+1,n+1); %Cure matrix initialized
end

clear count
for count=1:te/recordstep
    for ii=1:m+1
        jj=(ii-1)*(n+1)+1:(ii-1)*(n+1)+1+n;
        Tm(ii,:,count)=T(jj,count);  %Temperature matrices populated at each time step
        
        alpham(ii,:,count)=alpha(jj,count); %Cure matrices populated at each time step
    end
end