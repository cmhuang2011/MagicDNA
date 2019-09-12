function NeighborM=Neibor3(input)
A=input;

n=size(A,1);
m=size(A,2);
A=wextend(2,'zpd',A,1);
B=zeros(n,m);
for i=1:n
   for ii=1:m
       if A(i+1,ii+2)>0   %right element
        B(i,ii)=  1;
       end
       if A(i+2,ii+1)>0     %down element
        B(i,ii)= 1;
       end     
       if A(i+1,ii)>0     %left element
        B(i,ii)=  1;
       end
       if A(i,ii+1)>0   %top element
        B(i,ii)=  1;
       end  
      
    
   end   
end
B(find(input~=0))=0;
B=B>0;
NeighborM=B;