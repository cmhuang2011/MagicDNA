
function   [plotH,HeadTail]= DrawStapp3(StapleCell,GetHyperB,DrawFlag,ax,pStapleH)    %plot 2D frame
CC=[];
ZZ=[];
if DrawFlag==0
    return
end

axes(ax) ; title('Staple 2D') ; hold on ;  % very important to hold on ; before plotting anything, the axe need to be set hold on;
set (gca,'Ydir','reverse') ;
% ax.Tag
% return
for hh=1:length(StapleCell)
    clysM=StapleCell{hh};
    clys=clysM(:,1);
    CC=union(CC,clys);
    ZZ=union(ZZ,clysM(:,2));
end

% nA=max(CC);

nA=max(CC) + (length(GetHyperB.containBundle))*2;  % spacing bundles

MaxZZ=max(ZZ);
nn=1;
mm=nA;
XCordinate=zeros(mm*nn,2);
YCordinate=zeros(mm*nn,2);
xxc=zeros(mm*nn,2);
yyc=zeros(mm*nn,2); 
zzc=zeros(mm*nn,2);
L=MaxZZ+10;
m=0;
for i=1:nA
    m=m+1;
    XCordinate(m,:) = [0 ,L]; %horizontal lines(X-direction)
    YCordinate(m,:) = [5*i ,5* i];
end
m2=0;
for j=1:nA
    m2=m2+1;
    xxc(m2,:)=[-5,-5];
    yyc(m2,:)=[5*j ,5* j];
end
%                  if DrawFlag==1
% %                 FFN=figure('Color',[1 1 1]);
% %                 set(FFN,'name','2D Stapple Panel');
% %
%                   for d=2:ceil(nA/10)      %Solid horizontal ref lines
% %                   plot([0 L]',[50*(d-1),50*(d-1)],style3{:})
%                   end
%
%                  end
style = {'color', [0.7,0.7,0.7]  ,'LineStyle','-','linewidth' ,0.5};
xc=XCordinate;
yc=YCordinate;
SaveXY = zeros(500,2) ; nS=1;
if DrawFlag~=0
    
    XXc =[xc' ;nan*ones(1,size(xc,1) )]  ; XXc= reshape(XXc , numel(XXc),1) ;
    YYc =[yc' ;nan*ones(1,size(yc,1) )]  ; YYc= reshape(YYc , numel(YYc),1) ;
    SaveXY(nS:nS+size(XXc,1)-1 ,:) =[XXc,YYc] ; nS=nS+size(XXc,1) ;
   
%     plot(xc', yc', style{:})  ;       %  hold on;
    style2 = {'color', 'blue','LineStyle','-','linewidth' ,2.5};
    plot([0 0]',[5,5*nA],style2{:},'HitTest','off' )
    
    %   style3 = {'color', 'black','LineStyle','-','linewidth' ,1.5};
    
    style4 = {'color', [0.7,0.7,0.7],'LineStyle','-','linewidth' ,0.5};
    for d=1:floor(L/20)      %Solid vertical ref lines
        xxx=[20*(d),20*(d),nan ]' ;
        yyy=[5,nA*5,nan ]' ;
        SaveXY(nS:nS+size(xxx,1)-1 ,:) =[xxx,yyy] ; nS=nS+size(xxx,1) ;
        
%         plot([20*(d),20*(d)]',[5,nA*5],style4{:})
        
%         text( 20*(d)  , -1,int2str(d*20),'FontSize',8,'color', 'black','clipping','on');
    end
    axis off equal
    set(gca, 'visible', 'off')
end
SaveXY=SaveXY(1:nS-1,:) ;
plot(SaveXY(:,1),SaveXY(:,2),style4{:},'HitTest','off' ) ;
diS=length(StapleCell) ;
plotH= cell(diS,1) ;



HeadTail=cell(1,4) ;
%         HeadTail=gobjects(diS,3) ;
CollectHeadTail =zeros(diS ,4) ;
for di=1:1:diS
    CellMat=  StapleCell{di};
%     
%     IndsSecond = CellMat(:,1) >24 ;  % hard code for ploting staples separate by bundles   
%     CellMat(IndsSecond,1 ) =  CellMat(IndsSecond,1 )  +2 ;
    
    [~,IndC5] = ismember( CellMat(:,1) ,GetHyperB.RelateTable(:,5) ) ;
    IndB = GetHyperB.RelateTable(IndC5,1) ;
     IndB(GetHyperB.RelateTable(IndC5,2)==-1 )= length(GetHyperB.containBundle)+1 ;
    CellMat(:,1 ) =  CellMat(:,1 )  + (IndB-1)*2 ;
    
    
    while 1
        colorindex=randi([0 5],1,3)     ;
        colorindex(:)=colorindex(:)/5;
        if ismember(sum(colorindex),[1 2])
            break
        end
    end
    %                     colorindex= pStapleH{di}.Color;    % if use colors as
    %                     3D
    style5 = {'color',[colorindex],'LineStyle','-','linewidth' ,2};
    
    %           style5 = {'color','r','LineStyle','-','linewidth' ,1.5};
    m=size(CellMat,1)/2;
    xxc=zeros(6*m-4,1);
    yyc=zeros(6*m-4,1);
    mi=0;
    %             xxc(1)=CellMat(1,2);
    %             yyc(1)=5*CellMat(1,1);
    for tt=1:size(CellMat,1)-1
        fff=mod(tt,2);
        switch fff
            case 1    %horizontal points
                mi=mi+1;
                xxc(mi)=CellMat(tt,2);
                yyc(mi)=5*CellMat(tt,1)+1  ;
            case 0    %vertical points
                
                CC=[5*CellMat(tt,1) ,5*CellMat(tt+1,1)];
                delta=0.1*abs(CC(1)-CC(2));
                JJ= xor((CellMat(tt-1,2)>CellMat(tt,2)),(CC(1)>CC(2)));
                if JJ==true
                    amp=(1.1)^delta;
                else
                    amp=-(1.1)^delta;
                end
                dy=linspace( CC(1),CC(2),5);
                
                dx=amp*sin(pi/abs(CC(1)-CC(2)).*(dy-CC(1)) );
                %                         if CellMat(tt,1)>CellMat(tt-1,1)
                %                             dy=flip(dy)
                %                             dx=flip(dx)
                %                         end
                dx(isnan(dx))=0;
                Rangedx= max(dx)-min(dx) ; meandx = mean(dx);
                dx=dx-meandx; dx=dx/Rangedx*2;
                dx=dx-dx(1) ;
                xxc(mi+1:mi+5)=  linspace( CellMat(tt,2),CellMat(tt+1,2),5)+dx   ;
                %                                         xxc(mi+1:mi+5)=  CellMat(tt,2)+dx   ;
                yyc(mi+1:mi+5)=  dy+1;
                mi=mi+5;
        end
    end
    mi=mi+1;
    xxc(mi)=CellMat(end,2);
    yyc(mi)=5*CellMat(end,1)+1;
    %           plot([CellMat(tt,2), CellMat(tt+1,2)]',[5*CellMat(tt,1), 5*CellMat(tt+1,1)],style5{:}) ;
    plotH{di} =plot(xxc,yyc,style5{:} ) ;
    %           HeadTail(di,1) = text(  CellMat(1,2)  ,  5*CellMat(1,1)+2,int2str(di),'FontSize',8,'color', 'red','clipping','on');
    %           HeadTail(di,2) = scatter( CellMat(1,2)  ,  5*CellMat(1,1)+1,'square','filled','r');
    %           HeadTail(di,3) =  scatter( CellMat(end,2)  ,  5*CellMat(end,1)+1,'d','filled','b');
%     HeadTail{di,1} = text(  CellMat(1,2)  ,  5*CellMat(1,1)+2, strcat('\leftarrow', int2str(di)) ,'FontSize',8,'color', 'red','clipping','on');
%     HeadTail{di,2} = scatter( CellMat(1,2)  ,  5*CellMat(1,1)+1,'square','filled','r');
%     HeadTail{di,3} =  scatter( CellMat(end,2)  ,  5*CellMat(end,1)+1,'d','filled','b');
%     HeadTail{di,4} =   text(  CellMat(end,2)  ,  5*CellMat(end,1)+2, strcat('\leftarrow',int2str(di)),'FontSize',8,'color', 'blue','clipping','on');   %---
    %-----------improve graphical efficiency June 24 2019
    CollectHeadTail(di,1:2) = [CellMat(1,2)  ,  5*CellMat(1,1)+1 ] ;
    CollectHeadTail(di,3:4) = [CellMat(end,2)  ,  5*CellMat(end,1)+1 ] ;
    
    
end
str= num2str([1:diS]');
% str= strcat('\leftarrow', num2str([1:diS]')) ;

HeadTail{1,1} = text( CollectHeadTail(:,1)' ,CollectHeadTail(:,2)', str,'FontSize',8,'color', 'red','clipping','on','HitTest','off'  );

HeadTail{1,2} = scatter( CollectHeadTail(:,1) ,CollectHeadTail(:,2) ,'square','filled','r','HitTest','off' ); 
HeadTail{1,3} = scatter( CollectHeadTail(:,3) ,CollectHeadTail(:,4) ,'d','filled','b','HitTest','off' );
HeadTail{1,4} = text( CollectHeadTail(:,3)' ,CollectHeadTail(:,4)', str,'FontSize',8,'color', 'blue','clipping','on','HitTest','off' );   %---



end

