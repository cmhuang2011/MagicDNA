function StapBaseList = InsertAllBaseForStaple( StapListCell  )
%% Corner notaion to all base
%   Detailed explanation goes here
StapBaseList=cell(size(StapListCell)) ;
for k=1: length(StapBaseList)
   Stp = StapListCell{k} ;
   
    QQ=abs(diff(Stp(:,2))) ;
    TotalBases = sum(QQ) + sum(QQ~=0) ;
    
   BasesAllocate =  zeros( TotalBases, 2) ;
   N_corner = size(Stp,1 ) ; n=1;
   for j_sec= 1: round(N_corner/2)
       
       if Stp(2*j_sec,2) > Stp(2*j_sec-1 ,2)
           Arr= [  Stp(2*j_sec,1)* ones(1+abs(diff(Stp(2*j_sec-1:2*j_sec ,2))) ,1) , (Stp(2*j_sec-1 ,2):1: Stp(2*j_sec ,2))'  ];
       else
           Arr= [  Stp(2*j_sec,1)* ones(1+abs(diff(Stp(2*j_sec-1:2*j_sec ,2))) ,1) , (Stp(2*j_sec-1 ,2):-1: Stp(2*j_sec ,2))'  ];  
       end
       BasesAllocate(n:n+size(Arr,1)-1 ,:) = Arr ;
       n=n+ size(Arr,1);

   end
    StapBaseList{k}  = BasesAllocate;
        
end



end

