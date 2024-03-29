function  AssignCustomSQ
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% [ HClattice ] = findHClattice( 1 ,[40 40]) ;
% [X,Y] = meshgrid(0:2:80,0:2:80) ;   % For larger Xsec, 07222020
[X,Y] = meshgrid(0:2:40,0:2:40) ;

SQlattice.SQcenter= [X(:) ,Y(:) ] ;

fH=figure(236); fH.Units='normalized';fH.OuterPosition=[0 0 1 1];  clf;
h_scatter =scatter(SQlattice.SQcenter(:,1) ,SQlattice.SQcenter(:,2),'o','filled');
h_scatter.CData=ones(length(SQlattice.SQcenter(:,1)),1 ) *[0,0,1] ;
CylinderIndex =[1:length(SQlattice.SQcenter(:,2))]' ;
hold on ;
handletext= text(SQlattice.SQcenter(:,1)+0.2 ,SQlattice.SQcenter(:,2),num2str([1:length(SQlattice.SQcenter(:,2))]' ),'HitTest','off','Visible','off' );
axis equal ; grid off ;

ax=gca;
ax.Position=[0.2 ,0.1, 0.6,0.8] ;

h_scatter.ButtonDownFcn=@(src,evn)ClickScatter(src,evn,CylinderIndex) ;

title('Left Click to select cylinders(Red). Right Click to cancel. Enter to Update  ','FontSize',16)
% handletext(1).HitTest
fH.NumberTitle = 'off' ;
fH.Name = 'custom SQ GUI' ;

r=1 ; N=60 ;
AngsDeg = linspace(0,360,N) ;
dxdy = r*[cosd(AngsDeg)' ,sind(AngsDeg)'] ;
% XYforCircumference = zeros(( size(HClattice.SQcenter,1)+1)*N ,2)  ;
    
XYforCircumference= repmat([dxdy; 0,0] , size(SQlattice.SQcenter,1) ,1 ) + repelem(SQlattice.SQcenter, N+1,1)  ;
XYforCircumference( (N+1):(N+1): end ,:) = ones(size(SQlattice.SQcenter,1) ,1)*[NaN, NaN] ;
plot(XYforCircumference(:,1) ,XYforCircumference(:,2) ,'k' , 'HitTest','off' );
%---------------
A_Inds =zeros(size(SQlattice.SQcenter,1) ,1);
A_Inds(1:2:end) = 1; A_Inds= A_Inds==1 ;

A_Centers =  SQlattice.SQcenter(A_Inds ,:) ;
B_Centers =  SQlattice.SQcenter(~A_Inds ,:) ;
XY_circum_A =  repmat(dxdy , 1 ,size(A_Centers,1) ) + repelem( reshape(A_Centers' , 1,numel(A_Centers)) , N, 1)  ;
XY_circum_B =  repmat(dxdy , 1 ,size(B_Centers,1) ) + repelem( reshape(B_Centers' , 1,numel(B_Centers)) , N, 1)  ;

ptH_A = patch(XY_circum_A(:,1:2:size(XY_circum_A,2)) ,XY_circum_A(:,2:2:size(XY_circum_A,2)) ,[0.8,0.8,0] ,'FaceAlpha',0.2) ;
ptH_B = patch(XY_circum_B(:,1:2:size(XY_circum_B,2)) ,XY_circum_B(:,2:2:size(XY_circum_B,2)) ,[0.2,0.7,0.7] ,'FaceAlpha',0.2) ;
ptH_A.FaceColor = 'flat'; ptH_B.FaceColor = 'flat';
ptH_A.CDataMapping='direct' ; ptH_B.CDataMapping='direct' ;
ptH_A.FaceVertexCData = repmat([0.8,0.8,0] ,size(ptH_A.Faces,1),1  );
ptH_B.FaceVertexCData = repmat([0.2,0.7,0.7],size(ptH_B.Faces,1),1  );

ShowText = uicontrol(fH,'Style','text','Units','normalized','Position',[0.85 0.1 0.15 0.1],'String','[5''->3'',3''->5'']','HorizontalAlignment','left','FontSize',16,'Enable','Inactive' );
str={'xSQ1';'xSQ2';'xSQ3';'xSQ4';'xSQ5'} ;
MultipleXSec = uicontrol(fH,'Style','popup','Units','normalized','Position',[0.85 0.25 0.1 0.1],'String',str ,'FontSize',16 );



uistack(ptH_A,'bottom') ; uistack(ptH_B,'bottom') ;
h_scatter.ButtonDownFcn=@(src,evn)ClickScatter(src,evn,CylinderIndex,ptH_A,ptH_B,A_Inds,ShowText ) ;
MultipleXSec.Callback=@(src,evn)SelectxSQ(src,evn,CylinderIndex, h_scatter ,ptH_A,ptH_B,A_Inds , ShowText ) ;

fH.KeyPressFcn= @(src,evn) Decide(src,evn,CylinderIndex, h_scatter ,MultipleXSec ,ptH_A,ptH_B,A_Inds )  ;

axis ij ;% y axis(row direction in cadnano) going down
axis equal ; 
%----------------------
            ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
            hLg= legend(ForLegend,'Click me for instructions','Location','northwest' ) ; hLg.String={'\bf{Click me for instructions}'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
            ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_customSQ( src,evn,ax);
            hLg.Units='normalized'; hLg.AutoUpdate ='off'; hLg.Position=[0.0063 0.9728 0.1569 0.025];

            ShowText.ButtonDownFcn=@(src,evn)uisetfont(ShowText);

end

function ClickScatter(src,evn ,CylinderIndex ,ptH_A,ptH_B,A_Inds,ShowText)

% d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2 ;
% Ind = find(d==min(d)) ;
% if evn.Button==1
%     src.CData(Ind,:) =[1,0,0] ;
% elseif evn.Button==3
%     src.CData(Ind,:) =[0,0,1] ;
% end
d=  (src.XData-evn.IntersectionPoint(1)).^2 + (src.YData-evn.IntersectionPoint(2)).^2 ;
Ind = find(d==min(d)) ;
if evn.Button==1
    src.CData(Ind,:) =[1,0,0] ;
    if A_Inds(Ind) ==0
        LocalBInd =sum(A_Inds(1:Ind)==0  ) ;
        ptH_B.FaceVertexCData(LocalBInd, : ) = [1,0,0];
    else
        LocalAInd =sum(A_Inds(1:Ind)==1  ) ;
        ptH_A.FaceVertexCData(LocalAInd, : ) = [1,0,0];
    end
elseif evn.Button==3
    
    src.CData(Ind,:) =[0,0,1] ;
    if A_Inds(Ind) ==0
        LocalBInd =sum(A_Inds(1:Ind)==0  ) ;
        ptH_B.FaceVertexCData(LocalBInd, : ) = [0.2,0.7,0.7];
    else
        LocalAInd =sum(A_Inds(1:Ind)==1  ) ;
        ptH_A.FaceVertexCData(LocalAInd, : ) = [0.8,0.8,0];
    end
    
end
IsRedDot = ismember(src.CData, [1,0,0],'rows' ) ;
% fprintf('N of each cylinder group = %i %i \n' ,[sum(A_Inds(IsRedDot)) , sum(~A_Inds(IsRedDot))]) ;
str='Cyl. direction' ;
str=[str newline '[5''->3'',3''->5'']'] ; sp='       ';sp1=' ';
str=[str newline sp1 num2str(sum(A_Inds(IsRedDot))) sp  num2str(sum(~A_Inds(IsRedDot)))  ] ;


ShowText.String = str ;
end

function SelectxSQ(src,evn,CylinderIndex, h_scatter ,ptH_A,ptH_B,A_Inds,ShowText ,MultipleXSec )
src.FontSize = ShowText.FontSize;

load('CustumSQLattice.mat','CustumSQLattice') ;
ThisHCLattice = CustumSQLattice{src.Value} ;

h_scatter.CData(ThisHCLattice,:) = ones( length(ThisHCLattice),1)*[1,0,0]  ;
NotSelect = setdiff(1:size(h_scatter.CData,1) ,ThisHCLattice ) ;
h_scatter.CData(NotSelect,:) =ones( length(NotSelect),1)*[0,0,1] ;

 ptH_B.FaceVertexCData= ones(size(ptH_B.FaceVertexCData,1),1 )*[0.2 0.7 0.7] ;
 ptH_A.FaceVertexCData= ones(size(ptH_A.FaceVertexCData,1),1 )*[0.8 0.8 0] ;

for k=1:length(ThisHCLattice)
Ind = ThisHCLattice(k) ;
    if A_Inds(Ind) ==0
        LocalBInd =sum(A_Inds(1:Ind)==0  ) ;
        ptH_B.FaceVertexCData(LocalBInd, : ) = [1,0,0];
    else
        LocalAInd =sum(A_Inds(1:Ind)==1  ) ;
        ptH_A.FaceVertexCData(LocalAInd, : ) = [1,0,0];
    end
    
end
IsRedDot = ismember(h_scatter.CData, [1,0,0],'rows' ) ;
% fprintf('N of each cylinder group = %i %i \n' ,[sum(A_Inds(IsRedDot)) , sum(~A_Inds(IsRedDot))]) ;
str='Cyl. direction' ;
str=[str newline '[5''->3'',3''->5'']'] ; sp='       ';sp1=' ';
str=[str newline sp1 num2str(sum(A_Inds(IsRedDot))) sp  num2str(sum(~A_Inds(IsRedDot)))  ] ;
ShowText.String =str ;


end



function Decide(src,evn,CylinderIndex, h_scatter ,MultipleXSec ,ptH_A,ptH_B,A_Inds )

% evn

SelectxSec = MultipleXSec.Value ; 
load('CustumSQLattice.mat','CustumSQLattice') ;
% ThisSQLattice = CustumSQLattice{src.Value} ;

% CustumSQLattice=[] ;
if strcmp(evn.Key, 'return' )
   
    Inds = find(ismember( h_scatter.CData, [1,0,0],'rows')) ;
    CustumSQLattice{SelectxSec} =  Inds ;
    %-----decide selected cylinders make sense or not
%     [ SQlattice ] = findHClattice( 1 ,[40 40]) ;
%     [X,Y] = meshgrid(0:2:80,0:2:80) ;  % For larger Xsec, 07222020
      [X,Y] = meshgrid(0:2:40,0:2:40) ;
    SQlattice.SQcenter= [X(:) ,Y(:) ] ;

    
    
    XY=[SQlattice.SQcenter(Inds,1),SQlattice.SQcenter(Inds,2)] ;

    OneBundle=  BundleCylinderHC(1,[],50*ones(1,size(XY,1)),100*ones(1,size(XY,1)),XY) ;
    HyperBundleWithOnlyOne=hyperbundle(1,OneBundle);  % first bundle 
    trial=0;
    while trial<=100
    [Doable,CellPairList,CellUnpair]=checkpairableH(HyperBundleWithOnlyOne,3);
     if  (Doable==1   ||  trial>=800)
           break;
     end
     trial=trial+1;
    end     
    Doable;
    %-------------
    if Doable==0  ||  size(XY,1)~=length(OneBundle.Zbase1 )
    opts = struct('WindowStyle','modal',... 
              'Interpreter','tex');
    f = errordlg('\fontsize{18} Cylinders can''t be paired or lumped. Consider different cross-section. Data is not updated','File Error',opts) ;
    
    else
    opts = struct('WindowStyle','modal',... 
              'Interpreter','tex');

    f = msgbox('\fontsize{18} Data Update',opts) ;
    fprintf('Number of Cylinders = %i \n',length(Inds) )
    save('CustumSQLattice.mat','CustumSQLattice') ;         
%     fprintf('Update selected HC Cylinders : [%s] \n',num2str(Inds')) ;
        
    end
      
end
% evn.Key

if strcmp(evn.Key, 'c' )

    h_scatter.CData= ones(size( h_scatter.CData,1 ),1 )*[0,0,1]  ;
    ptH_B.FaceVertexCData = ones(size( ptH_B.FaceVertexCData,1 ),1 )*[0.2,0.7,0.7] ;
    ptH_A.FaceVertexCData = ones(size( ptH_A.FaceVertexCData,1 ) ,1) *[0.8,0.8,0];
%     SDF=3
end



end

