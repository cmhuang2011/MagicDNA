function OverhangInitial( src,evn )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

ax= findobj(gcf,'Tag','MechOH3D');
ax2= findobj(gcf,'Tag','MechOH2D');

%             GO_inthistab= findobj( findobj(0,'Tag','ss_Assembly'),'Type','UIcontrol' ) ;
%             for k=1:length(GO_inthistab)
%                 if ~isfield(GO_inthistab(k).UserData,'keepme');   delete(GO_inthistab(k)) ;   end
%             end

% ax= findobj(0,'Tag','AssemblyMain');
axes(ax); cltab ;
axes(ax2);  imshow([pwd filesep 'images' filesep 'OH.jpg']) ;axis off ;

axes(ax); axis equal;
ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;

coeff=1;
patchH=cell(length( GetHyperB.containBundle),1) ;

% scatH=cell(length( GetHyperB.containBundle),1) ;
% save nicks : [Bundle Cyl, Base1,Base2, Gx,Gy,Gz];
ExcludeEdge = 6;   % unit: nt, bases from two side, hard code
% ExcludeEdge = 11;   % unit: nt, bases from two side   MagicDNA original
% setting

AllNicks=[0,0,0,0,0,0,0];
GatherHelix=zeros(round(1.2*GetHyperB.EstimateDsScafL) ,3 ) ; cG=1;

SQpattern= [ 2 4 7 10 12 15 18 20 23 26 28 31 ] ; SQpattern=[7 15 23 31] ;
HCpattern= [2 6 9 13 16 20] ; HCpattern=[6 13 20] ;

for kk=1:length( GetHyperB.containBundle)
    QQ=GetHyperB.containBundle{kk}.HelixXYZGStapNoPM(1,1) ;
    %         QQ=GetHyperB.containBundle{kk}.HelixXYZGStapNoPM(1,1) ;
    
    QQScaf=GetHyperB.containBundle{kk}.HelixXYZGStapNoPM(1,0) ;
    for cyli=1:length(QQ)
        %         pp=plot3(QQ{cyli}(:,1)/coeff,QQ{cyli}(:,2)/coeff,QQ{cyli}(:,3)/coeff,'-k','LineWidth',3 ,'HitTest','off') ;
        %         pp=line(QQ{cyli}(:,1)/coeff,QQ{cyli}(:,2)/coeff,QQ{cyli}(:,3)/coeff,'LineWidth',3 ) ;
        %         pp.Color=[0.1,0.1,0.9,0.7];
        %         pp.Color=[0.8,0.8,0.8,0.7];
        GatherHelix( cG:cG+length(QQ{cyli}(:,1)) , :) = [[QQ{cyli}(:,1)/coeff,QQ{cyli}(:,2)/coeff,QQ{cyli}(:,3)/coeff]; nan nan nan] ; cG=cG+length(QQ{cyli}(:,1))+1;
    end
    allpoints =[cell2mat(QQ);cell2mat(QQScaf)] ;
    shp=alphaShape(allpoints,1.5);   %2
    patchH{kk}=plot(shp);   %plot Volume
    TargetFaces= 1200;
    R=TargetFaces/size(patchH{kk}.Faces ,1);
    %                 reducepatch(patchH{kk},R) ;
    patchH{kk}.EdgeColor='none' ;
    patchH{kk}.FaceAlpha=0.1;
    patchH{kk}.PickableParts ='none';
    patchH{kk}.FaceColor= [0,0.2,0.8 ];
    set( patchH{kk},'HitTest','off') ;
    set( patchH{kk},'PickableParts','none')    ;
    
    %     colR=[2*rand(1,1),2*rand(1,1),rand(1,1)] ;
    %     colR=colR/max(colR) ;
    %     patchH{kk}.FaceColor=colR ;
    
    QQ2=GetHyperB.containBundle{kk}.HelixXYZGStapNoPM(1.5,1) ;
    %         QQ2=GetHyperB.containBundle{kk}.HelixXYZGStapNoPM(1.5,1) ;
    if strcmp( GetHyperB.containBundle{kk}.type, 'HC')
        Inv=1 ; Per=21;  ppatern= HCpattern ;
    else
        Inv=1;  Per=32;   ppatern= SQpattern ;     % density of staple overhang position, fixing indexing issue 05202019
    end
    %     QQ2=QQ2(10:end-10,:) ;
    %     scatH{kk} =cell( )
    for cyli=1:length(QQ2)
        %         if ( kk==1 && ismember(cyli,[1:3:18]) ) ||  ( kk==2 && ismember(cyli,[3:3:18]) )  % hard code
        
        QQ2{cyli}=QQ2{cyli}(ExcludeEdge:end-ExcludeEdge,:) ;
        Inds=GetHyperB.containBundle{kk}.Zbase1(cyli)+ExcludeEdge-1:GetHyperB.containBundle{kk}.Zbase2(cyli)-ExcludeEdge;
        Inds=Inds(1:Inv:end-1) ;
        QQ2{cyli}=  0.5*(QQ2{cyli}(1:Inv:end-1,:)+  QQ2{cyli}(2:Inv:end,:)) ;
        %         plot3(QQ2{cyli}(:,1)/coeff,QQ2{cyli}(:,2)/coeff,QQ2{cyli}(:,3)/coeff,'--','LineWidth',1)
        InorOut= inShape(shp,QQ2{cyli}(:,1),QQ2{cyli}(:,2),QQ2{cyli}(:,3));   %1-> in, 0->out
        %          scatter3(QQ2{cyli}(~InorOut,1),QQ2{cyli}(~InorOut,2),QQ2{cyli}(~InorOut,3),'o'  ) ;
        
        
        %          BasesL=Inds(~InorOut)' ;
        IndPerpendicular= and(~InorOut,  ismember(mod(Inds',Per),ppatern) ) ;
        
        
        BasesL=Inds(IndPerpendicular)' ;
        BasesH=BasesL+ones(size(BasesL));
        
        Saves=[kk*ones(size(BasesL)),cyli*ones(size(BasesL)),BasesL,BasesH,QQ2{cyli}(IndPerpendicular ,:)  ] ;
        %          AllNicks
        if ~isempty(Saves)
            AllNicks=union(AllNicks,Saves,'rows') ;
        end
        %          plot3(QQ2{cyli}(:,1),QQ2{cyli}(:,2),QQ2{cyli}(:,3),'.-r')
        %           scatter3(QQ2{cyli}(8,1),QQ2{cyli}(8,2),QQ2{cyli}(8,3),'ok','filled')
        %           scatter3(QQ2{cyli}(9,1),QQ2{cyli}(9,2),QQ2{cyli}(9,3),'ob','filled')
        %
        %    [cyli sum(IndPerpendicular)]
    end
end

pp=line(GatherHelix(:,1), GatherHelix(:,2) ,GatherHelix(:,3),'LineWidth',2) ;  pp.Color=[0.8,0.5,0.2,0.7];

AllNicks=setdiff(AllNicks,zeros(1,7),'rows');

ax.UserData.NicksScatH=  scatter3(AllNicks(:,5),AllNicks(:,6),AllNicks(:,7),56 ,'ok' ,'filled') ;
ax.UserData.NicksScatH.CData=zeros(size(AllNicks,1),3);
ax.UserData.NicksScatH.ButtonDownFcn=@(src,evn)showNicks(src,evn,AllNicks)  ;

ax.UserData.AllNicks=AllNicks;

btnAdd=findobj(gcf,'Tag','OHAdd'); btnAdd.Enable='on' ;
btnAdd3=findobj(gcf,'Tag','ApplyAll'); btnAdd3.Enable='on' ;

btnAdd.Callback= @(src,evn)addOH(src,evn) ;
btnAdd3.Callback= @(src,evn)ApplyAll(src,evn) ;

if ~strcmp(src.String,'Clear')
    ax.UserData.Conn=cell(200,1) ;
    ax.UserData.text=cell(200,2) ;
end
src.String='Clear';
t=findobj(gcf,'Tag','OHTable');
t.Data=''; t.UserData='';
t.CellEditCallback=@(src,evn)UpdateConn(src,evn) ;
% t.Position(2)=0.25;
% t.CellSelectionCallback=@(src,evn)CellSelection(src,evn) ;
% t.ButtonDownFcn =@(src,evn)btDown(src,evn) ;
h1 = light('Position',[-100 100 100],'Style','infinite'); h2 = light('Position',[20 -20 20],'Style','infinite');

xlabel('X') ; ylabel('Y') ; zlabel('Z') ;
end

function ApplyAll(src,evn)
t=findobj(gcf,'Tag','OHTable');
if ~isempty(t.Data)
    if size(t.Data ,1 )>1
        for k=2: size(t.Data ,1 )
            t.Data{k,2}  = t.Data{1,2} ;
            t.Data{k,3}  = t.Data{1,3} ;
            t.Data{k,5}  = t.Data{1,5} ;
            t.Data{k,6}  = t.Data{1,6} ;
            t.Data{k,7}  = t.Data{1,7} ;
            t.Data{k,8}  = t.Data{1,8} ;
            t.Data{k,9}  = t.Data{1,9} ;
        end
        eee.Indices=[1,3 ] ;
        UpdateConn(t,eee)  ;
        eee.Indices=[1,9 ] ;
        UpdateConn(t,eee)  ;
        
    end
end
end


function  showNicks(src,evn,AllNicks)
evn.IntersectionPoint;

XYZ(:,3)= src.ZData;
XYZ(:,1)= src.XData;
XYZ(:,2)= src.YData;
distances = sqrt(sum(bsxfun(@minus, XYZ, evn.IntersectionPoint).^2,2));
ind=find(distances==min(distances));

% scatter3(XYZ(ind,1),XYZ(ind,1),XYZ(ind,3)
% if ~ismember([1,0,0],src.CData,'rows') &&  evn.Button==1
% src.CData(ind,:) = [1,0,0];
% fprintf('AllNicksRR = %d %d %d %d %d %d %d \n',AllNicks(ind,:)) ;
% src.UserData.red=AllNicks(ind,:) ;
% elseif  ~ismember([0,0,1],src.CData,'rows') &&  evn.Button==3
% src.CData(ind,:) = [0,0,1];
% fprintf('AllNicksBB = %d %d %d %d %d %d %d \n',AllNicks(ind,:)) ;
% src.UserData.blue=AllNicks(ind,:) ;
%
% end
OHtext1=findobj(src.Parent.Parent,'Tag','OHText1') ; OHtext1.FontSize=12 ; OHtext1.Position(3)=0.14 ;
OHtext2=findobj(src.Parent.Parent,'Tag','OHText2') ; OHtext2.FontSize=12 ;  OHtext2.Position(3)=0.14 ;

if   evn.Button==1
    [Exist,ind2]=ismember([1,0,0],src.CData,'rows');
    if Exist==1;src.CData(ind2,:) = [0,0,0];end
    if sum(src.CData(ind,:) ==[0.5,0,0]) ~=3  && sum(src.CData(ind,:) ==[0,0,0.5]) ~=3
        src.CData(ind,:) = [1,0,0];
        
        fprintf('Selected Red dot = %d %d %d %d  \n',AllNicks(ind,1:4)) ;
        SStr= sprintf('Selected Red  = %d %d %d %d',AllNicks(ind,1:4))  ;
        OHtext1.String = SStr;
        src.UserData.red=AllNicks(ind,:) ;
    end
elseif   evn.Button==3
    [Exist,ind2]=ismember([0,0,1],src.CData,'rows');
    if Exist==1;src.CData(ind2,:) = [0,0,0];end
    if sum(src.CData(ind,:) ==[0.5,0,0]) ~=3  && sum(src.CData(ind,:) ==[0,0,0.5]) ~=3
        src.CData(ind,:) = [0,0,1];
        fprintf('Selected Blue dot = %d %d %d %d  \n',AllNicks(ind,1:4)) ;
        SStr2= sprintf('Selected Blue  = %d %d %d %d',AllNicks(ind,1:4))  ;
        OHtext2.String = SStr2;
        
        src.UserData.blue=AllNicks(ind,:) ;
    end
end
% axis equal;
end

function addOH(src,evn)
t=findobj(gcf,'Tag','OHTable');
ax=gca ;
BAori=     [ones(1,3); 0.94*ones(1,3)] ;
%     0.9400    0.9400    0.9400
%check blue and red dot exist
if ismember([1,0,0], ax.UserData.NicksScatH.CData,'rows') && ismember([0,0,1], ax.UserData.NicksScatH.CData,'rows')
    [~,RedInd] = ismember([1,0,0], ax.UserData.NicksScatH.CData,'rows')  ;
    [~,BlueInd] = ismember([0,0,1], ax.UserData.NicksScatH.CData,'rows')  ;
    
    t.Visible='on'   ;
    % sdfsf=3
    n=size(t.Data,1) ;
    Redstr= num2str(ax.UserData.NicksScatH.UserData.red(1:4)) ;  RedLoc=ax.UserData.NicksScatH.UserData.red(5:7) ;
    Bluestr= num2str(ax.UserData.NicksScatH.UserData.blue(1:4)) ;  BlueLoc=ax.UserData.NicksScatH.UserData.blue(5:7) ;
    
    if n==0
        t.Data={'A1','3''',true,'B1','3''',8,8,'Connected',true};   %initialize
        t.UserData={Redstr,'3''',true,Bluestr,'3''',[RedInd, BlueInd] };
        %     ax.UserData.text{1,1}=text(RedLoc(1),RedLoc(2),RedLoc(3),'\leftarrow A1' ,'Color','red','FontSize',14) ;
        %     ax.UserData.text{1,2}=text(BlueLoc(1),BlueLoc(2),BlueLoc(3),'\leftarrow B1','Color','blue','FontSize',14) ;
        MM=[RedLoc;BlueLoc];
        ax.UserData.Conn{1} = plot3(MM(:,1),MM(:,2),MM(:,3),'LineWidth',2) ;
        ax.UserData.Conn{1}.ButtonDownFcn=@(src,evn)CancelOH(src,evn,ax) ;
        t.BackgroundColor=BAori(1,:) ;
    else
        CellData=    t.Data; CurrenN= str2num(t.Data{n,1}(2:end)) ;
        CellData{n+1,1}=strcat('A',num2str(CurrenN+1));    CellData{n+1,2}='3''';    CellData{n+1,3}=true;
        CellData{n+1,4}=strcat('B',num2str(CurrenN+1));    CellData{n+1,5}='3''';
        CellData{n+1,6}=8;     CellData{n+1,7}=8;     CellData{n+1,8}='Connected';  CellData{n+1,9}=true;
        
        
        t.Data=CellData;
        
        %     ax.UserData.text{n+1,1}=text(RedLoc(1),RedLoc(2),RedLoc(3),strcat('\leftarrow' ,'A',num2str(n+1)) ,'Color','red' ,'FontSize',14) ;
        %     ax.UserData.text{n+1,2}=text(BlueLoc(1),BlueLoc(2),BlueLoc(3),strcat('\leftarrow' ,'B',num2str(n+1)),'Color','blue','FontSize',14) ;
        %     ax.UserData.text{1,1}.Position=[ ax.UserData.text{1,1}.Position ;[RedLoc(1),RedLoc(2),RedLoc(3)] ];
        %     ax.UserData.text{1,1}.String={ ax.UserData.text{1,1}.String ;  strcat('\leftarrow' ,'A',num2str(n+1)) };
        
        MM=[RedLoc;BlueLoc];
        ax.UserData.Conn{n+1} = plot3(MM(:,1),MM(:,2),MM(:,3),'LineWidth',2) ;   % save plot handle here.
        ax.UserData.Conn{n+1}.ButtonDownFcn=@(src,evn)CancelOH(src,evn,ax) ;
        CellData=    t.UserData;
        CellData{n+1,1}=Redstr;    CellData{n+1,2}='3''';    CellData{n+1,3}=true;
        CellData{n+1,4}=Bluestr;    CellData{n+1,5}='3'''; CellData{n+1,6}=[RedInd, BlueInd] ;
        t.UserData=CellData;
        
        BaColor2 =  repmat(BAori,n,1); BaColor2=BaColor2(1:n+1,:);
        %     BaColor=BaColor(1:n+1,:) ;
        % %     t.UserData.oriBaColor = BaColor ;
        t.BackgroundColor=BaColor2 ;
        
    end
    ax.UserData.NicksScatH.CData(RedInd,:) = [0.5,0,0] ;
    ax.UserData.NicksScatH.CData(BlueInd,:) = [0,0,0.5] ;
    
    %     ax.UserData.NicksScatH.CData=zeros(size(ax.UserData.NicksScatH.CData))
    %     ;  %
else
    
    warndlg('Use right and left clicks to select overhangs locations','Warning');
end

end

function CancelOH(src,evn,ax)
t=findobj(gcf,'Tag','OHTable');
BAori=     [ones(1,3); 0.94*ones(1,3)] ;
IndsValid= cellfun(@isempty,ax.UserData.Conn); IndsValid=find(~IndsValid );
for k=1:length(IndsValid)
    if isequal(src,ax.UserData.Conn{k} )
        NthConn=k ;
    end
end
if  src.LineWidth ==2 % hightlight to cancel
    %    for
    if evn.Button==1
        
        opts.Interpreter = 'tex';
        % Include the desired Default answer
        opts.Default = 'Yes';
        % Use the TeX interpreter to format the question
        quest = 'Are you sure to cancel this ??';
        answer = questdlg(quest,'Cancel this pair of overhangs',...
            'Yes','No',opts) ;
        
        if strcmp(answer,'Yes')
            ax.UserData.NicksScatH.CData(t.UserData{NthConn,6} , : )   =zeros(2,3) ;
            delete( ax.UserData.Conn{NthConn}) ;
            ax.UserData.Conn(NthConn)=[];
            t.UserData(NthConn,:) =[];
            t.Data(NthConn,:) =[];
            t.BackgroundColor(NthConn,:)=[];
        end
    else
        if evn.Button==2
            src.LineWidth =0.3 ;% change to thin
            t.Data{NthConn,9}=false;
        end
        
    end
else
    %         evn
    if evn.Button==1
        BaColor2 =  repmat(BAori,k,1); BaColor2=BaColor2(1:k+1,:);
        
        t.BackgroundColor=BaColor2;
        t.BackgroundColor(NthConn,:)=[0.9,0.9,0.6] ;
    elseif evn.Button==3
        src.LineWidth =2 ;% change to bold, ready to be able to delete.
        t.Data{NthConn,9}=true;
        
    end
    
end

end


function UpdateConn(src,evn)
if  evn.Indices(2)==3  % only for chaning enabability
    ax= findobj(gcf,'Tag','MechOH3D');
    Conn=ax.UserData.Conn;
    for kk=1: size(src.Data,1)
        if src.Data{kk,3}==0
            Conn{kk}.LineStyle=':' ;
        else
            Conn{kk}.LineStyle='-' ;
        end
    end
    %     axis equal;
end

if  evn.Indices(2)==9  % only for chaning enabability
    ax= findobj(gcf,'Tag','MechOH3D');
    Conn=ax.UserData.Conn;
    for kk=1: size(src.Data,1)
        if src.Data{kk,9}==0
            Conn{kk}.LineWidth=0.3 ;
        else
            Conn{kk}.LineWidth=2 ;
        end
    end
    %     axis equal;
end


end


function CellSelection(src,evn)
% evn.Indices;
% UpdateConn(src,evn) ;
% src.ColumnFormat{8}
% if ~isempty(evn.Indices)
%     if evn.Indices(2) ==8
%      RowInd= evn.Indices(1) ;
%         if strcmp( src.Data{RowInd,2} ,src.Data{RowInd,5})
%              src.Data{RowInd,8} ='Ar_Be' ;
%             src.ColumnFormat{8} = {'Ar_Be','Ae_Br'} ;
%
%
%         else
%                src.Data{RowInd,8} ='r_r' ;
%             src.ColumnFormat{8} = {'r_r','e_e'} ;
%
%         end
%         drawnow;
%
%
%     end
% else
%
% %    fprintf('empty CellSelection \n ') ;
%
% end
fprintf(' CellSelection  %i\n ',randi(10)) ;
drawnow;
% src.ColumnFormat{8}

end


function btDown(src,evn)

fprintf(' btDown \n ') ;

end





