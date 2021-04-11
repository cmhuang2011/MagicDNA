function [ TTAllTraj ] = FcnCancelDrift( AllTraj,boxs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


%% cancel drfit out of box for trajectory file ;

boxsize= boxs ;

boxCenter = 0.5*boxsize*ones(1,3) ;
% figure; 
% hx=histogram(AllTraj(:,1,:)); hold on ;
% hy=histogram(AllTraj(:,2,:)); hz= histogram(AllTraj(:,3,:)) ; 
% 
% legend([hx,hy,hz], 'x ' ,'y ' ,'z '  ) ;
% tic
for iF=1: size(AllTraj,3)
    ConfT= AllTraj(:,1:3,iF) ;
 Center = mean(ConfT(:,1:3)) ;
 ConfT=ConfT - repmat(Center, size(ConfT,1),1) + repmat(boxCenter, size(ConfT,1),1) ;
 AllTraj(:,1:3,iF) =ConfT ;
end


TTAllTraj=AllTraj ;
for iF=192: size(AllTraj,3)
    ConfT= TTAllTraj(:,1:3,iF)- repmat(0.5*ones(1,3)*boxs, size(TTAllTraj,1) ,1)    ;
    
%       figure(7);  clf;
% hx1=histogram(ConfT(:,1)); hold on ;
% hy1=histogram(ConfT(:,2)); hz1= histogram(ConfT(:,3)) ; 
%   legend([hx1,hy1,hz1], 'x ' ,'y ' ,'z '  ) ;
%    eedge= 0:1:boxsize ;
   eedge= -1*boxsize:5:2*boxsize ;

[N1,edges1] = histcounts(ConfT(:,1) ,eedge ) ;
[N2,edges2] = histcounts(ConfT(:,2) ,eedge) ;
[N3,edges3] = histcounts(ConfT(:,3)  ,eedge) ;


  % x 
   [L,num1] = bwlabel(N1) ;
   if  num1==2
       
    IndMaxInFirst = find(L==1,1,'last') ;
    ThresHold= edges1(IndMaxInFirst+2) ;

    IndFs=ConfT(:,1)>ThresHold ;
    ConfT(IndFs,1)=  ConfT(IndFs,1)-boxsize;    
    meanx=mean(ConfT(IndFs,1)) ;
    ConfT(:,1)=   ConfT(:,1) -meanx +boxsize/2 ;

   end
   %
  % y 
   [L,num2] = bwlabel(N2) ;
   if  num2==2
       
    IndMaxInFirst = find(L==1,1,'last') ;
    ThresHold= edges2(IndMaxInFirst+2) ;

    IndFs=ConfT(:,2)>ThresHold ;
    ConfT(IndFs,2)=  ConfT(IndFs,2)-boxsize;    
    meany=mean(ConfT(IndFs,2)) ;
    ConfT(:,2)=   ConfT(:,2) -meany +boxsize/2 ;

   end
   %   
   % z 
   [L,num3] = bwlabel(N3) ;
   if  num3==2

    IndMaxInFirst = find(L==1,1,'last') ;
    %       ThresHold= edges3(IndMaxInFirst+1) ;
    ThresHold= mean( edges3([IndMaxInFirst+1,IndMaxInFirst+2])) ;

    IndFs=ConfT(:,3)>ThresHold ;
    ConfT(IndFs,3)=  ConfT(IndFs,3)-boxsize;  
    meanz=mean(ConfT(IndFs,3)) ;
    ConfT(:,3)=   ConfT(:,3) -meanz +boxsize/2 ;


    [Nz,edgesZ] = histcounts(ConfT(:,3)) ;
    [Lz,numz] = bwlabel(Nz) ;

%       if numz==2
%         
%          figure(38);  subplot(2,1,1);
%          histogram( AllTraj(:,3,iF)); 
%         subplot(2,1,2);
%           histogram( ConfT(:,3))
%           sdfsf=3 
%          
%       end
      
       
   end
   %     
   
   TTAllTraj(:,1:3,iF)=ConfT ;
    
%    if iF==161
%     figure(72);  clf ;
%     hx1=histogram(ConfT(:,1)); hold on ;
%     hy1=histogram(ConfT(:,2)); hz1= histogram(ConfT(:,3)) ; 
%     legend([hx1,hy1,hz1], 'x ' ,'y ' ,'z '  ) ;
%    end
%    
   
end

% toc









end

