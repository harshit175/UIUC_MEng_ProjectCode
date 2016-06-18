function [CT ] = modifymatrices( CT,b )
%MODIFYMATRICES Summary of this function goes here
%Making changes to CT, DT and K matrix values

%the code doesn't change the size of the stiffness matrix, rather rows and 
% columns corresponding to constant temperature are made zero. 
%diagonal elements in rows corresponding to constant temperature are made 1


%rows and columns made zero
CT(:,b)=0;
CT(b,:)=0;

%Diagonal elements made 1
for q=1:length(b)
    
    CT(b(q),b(q))=1; 
    
end
end

