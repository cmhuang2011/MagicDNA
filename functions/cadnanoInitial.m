function cadnanoInitial( src,evn,jsonSlider1,jsonSlider2,jsonPop ,jsonPop2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% clc
% profile on
% profile on  % last optimized, 08/26/2019

tic;
fH=gcf;

fprintf('initial cadnano \n') ;

ax= findobj(fH,'Tag','json3D') ;
ax2=findobj(fH,'Tag','json2D') ;
btn2= findobj(fH,'Tag','btn_json2') ;
btn3= findobj(fH,'Tag','btn_json3') ;
btn4= findobj(fH,'Tag','btn_json4') ; btn4.Enable='off' ;

JsonTxt_scaf =findobj(gcf, 'Tag','JsonTxt_scaf') ;
JsonTxt_stap =findobj(gcf, 'Tag','JsonTxt_stap') ;

% JsonTxt_scaf.ButtonDownFcn = @(src,evn)JsonTxt_scafButtonDown(src,evn) ;
t_json=findobj(fH,'Tag','t_json') ;

axes(ax) ;
cltab; hold on;
axes(ax) ;
ss_Assembly= findobj(fH,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;
skipBase= GetHyperB.skipBase ;
insertBase= GetHyperB.insertBase ;
GetHyperB.ObjFigure = fH ;

axes(ax);
hold on; axis  equal;


%----------2D
% GetHyperB.drawBCBdata({GetHyperB.ScafRouting},ax2);      %scaffold
axes(ax2) ; ax2.UserData.NoZ=1 ;
[plotH,HeadTail]=DrawStapp3(GetHyperB.StapList3,GetHyperB,1,ax2,[]) ;            % 2D staple
% HeadTail: handle of [text] [scatter] [scatter] [text]

answer2 = questdlg('How to color the staples?', ...
    'Staple Color', ...
    'By 5'' in bundles','By Staple graph','By 5'' in bundles');
UseStapleGraph = 0;
switch answer2
    case 'By 5'' in bundles'
        TurbulanceColor= 0.1*rand(length(GetHyperB.containBundle),3) + 0*ones(length(GetHyperB.containBundle),3 ) ;  %
        CC=  get(groot,'defaultAxesColorOrder') ; CC(1,:)=[];
        CC =repmat(CC ,ceil(length(GetHyperB.containBundle)/size(CC,1) ),1 ) ;
        CC=CC(1:length(GetHyperB.containBundle) ,: );
        BundleColorForStaples=CC+TurbulanceColor ; BundleColorForStaples(BundleColorForStaples>1)=1;
        
        for k=1:length(GetHyperB.StapList3 )
            HeadC5Cylinder =GetHyperB.StapList3{k}(1,1) ;
            BundleInd = unique(GetHyperB.RelateTable( GetHyperB.RelateTable(:,5)==HeadC5Cylinder ,1   ) ) ;
            plotH{k}.Color = BundleColorForStaples(BundleInd,: ) ;  % bundle's staples color
            plotH{k}.HitTest='off' ;
        end
        h_bingraph=[];
   
        
    case 'By Staple graph'
        UseStapleGraph=1 ;
        [ h_bingraph,Graph,ST_list ]= GetHyperB.FindStapGraph ;
        axes(ax2);
        bins =  conncomp(Graph,'OutputForm','vector','Type','weak') ;
        %         [edgebins,iC] = conncomp(Graph,'Type','weak')
        edgebinsWeak = conncomp(Graph,'OutputForm','cell','Type','weak') ;
        ST =ST_list;
        [table2array(Graph.Edges) , ST_list]; SSSTTT=table2array(Graph.Edges) ;
        EdgeGroups = zeros(size(ST,1 ) ,1) ;
        for k = 1 : length(GetHyperB.StapList3 )
            [EdgeInd,~] = find(ST==k) ;
            if isempty(EdgeInd)
            else
                %                 length(EdgeInd)>1
                EdgeIndOne= EdgeInd(1) ;
                Onenode =  ST(EdgeIndOne,1) ;
                EdgeGroups(EdgeIndOne) =  bins(Onenode) ;
                %             else
            end
        end
        h_bingraph.EdgeCData =EdgeGroups;
        h_bingraph.LineWidth =2;
        ST2=[ST ,  bins(ST(:,1))', bins(ST(:,2))'  , bins(SSSTTT(:,1))', bins(SSSTTT(:,2))' ] ;
        h_bingraph.UserData.Graph=   Graph ;
        
        
        TurbulanceColor= 0.2*rand(max(bins),3) + 0*ones(max(bins),3 ) ;  %
        CC=  get(groot,'defaultAxesColorOrder') ; CC(1,:)=[];
        CC =repmat(CC ,ceil(max(bins)/size(CC,1) ),1 ) ;
        CC=CC(1:max(bins) ,: );
        BinsColorForStaples=CC+TurbulanceColor ;
        BinsColorForStaples(BinsColorForStaples>1)=1;
        
        EdgesColors = zeros(length(EdgeGroups),3) ;
        for k=1: max(bins)  % index of Group
            NodesInThisGroup = edgebinsWeak{k} ;
            EdgeList = table2array(Graph.Edges) ;
            
            
            IndEdge1 = ismember( EdgeList(:,1),NodesInThisGroup) ; %find(IndEdge1)
            
            EdgesColors(IndEdge1,:) = ones(sum(IndEdge1),1) *BinsColorForStaples(k,:);
        end
        
        h_bingraph.EdgeColor =EdgesColors;
        for bi = 1:max(bins)
            Staps =  find(bins==bi) ;
            if length(Staps)==1
                plotH{Staps}.Color =[0.8,0.8,0.8] ;  % single staple color
                plotH{Staps}.HitTest='off' ;
            else
                for k=1: length(Staps)
                    plotH{Staps(k)}.Color =BinsColorForStaples(bi,:) ;  % single staple color
                    plotH{Staps(k)}.HitTest='off' ;
                end
            end
        end
        h_bingraph.ButtonDownFcn=@(src,evn)StapGraphButtonDown(src,evn) ;
end
% return


% ax2=findobj(0,'Tag','json2D')
for k=1:length(plotH); plotH{k}.HitTest='off' ; end
% for k=1:length(HeadTail); HeadTail{k}.HitTest='off' ; end
% for k=1:length(HeadTail); HeadTail(k).HitTest='off' ; end


% return

%------2D scaffold
surfH= cell( length( GetHyperB.scafC5), 1)  ;  % scaffold 2D surf graphic handles in cell
ScafForScatterAll = []; ScafForScatterAllSpacing= [] ;
ScafForScatterIndividual = cell( length( GetHyperB.scafC5), 1)  ;  %scaffold 2D surf graphic handles in cell, later inserting spacing ~= C5

for scaf_j = 1 : length( GetHyperB.scafC5)
    CellMat= GetHyperB.scafC5{scaf_j} ;
    BaseRoute = interpolateBase( CellMat ) ;
    % CellMat = interpolateBase( CellMat ) ;   % test
    
    ScafForScatter= setdiff(BaseRoute , skipBase,'rows' ,'stable') ;  %C5 notation
    %-------introduce insert
    ScafForScatterExtraCol = [ScafForScatter , [1: size(ScafForScatter,1)]']; 
    [~,bb] = ismember(insertBase,ScafForScatter ,'rows'   ) ;
    
    insertBaseExtraCol = [insertBase , bb+0.5]; 
    Merge= [ScafForScatterExtraCol; insertBaseExtraCol] ;
    sortMerge = sortrows(Merge,3) ;
    ScafForScatter=sortMerge(:,1:2) ;
    
%     insertBase;
%     sdfsf=3
    %--------
    if scaf_j==1
        ScafForScatterAll=   ScafForScatter ;
    else
        ScafForScatterAll=   [ScafForScatterAll; ScafForScatter] ;
    end
    
    
    ScafForScatterOri=ScafForScatter ;  % save 2D scaffould base routing without spacing
    [~,IndC5] = ismember( ScafForScatter(:,1) ,GetHyperB.RelateTable(:,5) ) ;
    IndB = GetHyperB.RelateTable(IndC5,1) ;
    IndB(GetHyperB.RelateTable(IndC5,2)==-1 )= length(GetHyperB.containBundle)+1 ;
    ScafForScatterIndividual{scaf_j} = ScafForScatter ;
    ScafForScatter(:,1 ) =  ScafForScatter(:,1 )  + (IndB-1)*2 ;  % insert spacing , here insert spacing so ~= C5
    if scaf_j==1
        ScafForScatterAllSpacing=   ScafForScatter ;
    else
        ScafForScatterAllSpacing=   [ScafForScatterAllSpacing; ScafForScatter] ;
    end
    
    [~,IndC5] = ismember( CellMat(:,1) ,GetHyperB.RelateTable(:,5) ) ;
    IndB = GetHyperB.RelateTable(IndC5,1) ;
    IndB(GetHyperB.RelateTable(IndC5,2)==-1 )= length(GetHyperB.containBundle)+1 ;
    CellMat(:,1 ) =  CellMat(:,1 )  + (IndB-1)*2 ;
    
    m=size(CellMat,1)/2;
    xxc=zeros(6*m-4,1);
    yyc=zeros(6*m-4,1);
    mi=0;
    %             xxc(1)=CellMat(1,2);
    %             yyc(1)=5*CellMat(1,1);
    for tt=1:size(CellMat,1)-1
        %      fff=mod(tt,2);
        if CellMat(tt,1)== CellMat(tt+1,1); fcase=1 ;
        elseif  CellMat(tt,2)== CellMat(tt+1,2); fcase=2 ;
        else;   fcase=3 ;
        end
        
        
        switch fcase
            case 1    %horizontal points
                mi=mi+1;
                xxc(mi)=CellMat(tt,2);
                yyc(mi)=5*CellMat(tt,1)-1  ;
            case 2    %vertical points
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
                
                xxc(mi+1:mi+5)=  CellMat(tt,2)-dx   ;
                yyc(mi+1:mi+5)=  dy-1;
                mi=mi+5;
            case 3
                CC=[5*CellMat(tt,1) ,5*CellMat(tt+1,1)];
                dy=linspace( CC(1),CC(2),5);
                xxc(mi+1:mi+5)= linspace( CellMat(tt,2),CellMat(tt+1,2),5);
                yyc(mi+1:mi+5)=  dy-1;
                mi=mi+5;
                
        end
    end
    %     toc
    mi=mi+1; %CellMat
    xxc(mi)=CellMat(end,2);
    yyc(mi)=5*CellMat(end,1)-1;
    xxc=xxc(1:mi) ; yyc=yyc(1:mi) ;
    % % % % plot(xxc,yyc,'-*' )
    x=xxc' ; y=yyc' ; z =zeros(size(y)) ;
    % % % x= CellMat(:,2)' ; y= 5*CellMat(:,1)' -1 ;  z =zeros(size(y)) ;
    
    % x=
    %     col = (1:length(xxc))*1000;  % This is the color
    col = ((1:length(xxc))-1)/( length(xxc)-1);  % This is the color
    
    surfH{scaf_j}=surface([x;x],[y;y],[z;z],[col;col] ,'facecol','no', 'edgecol','interp','LineWidth',1.5    );
    % surfH=surface([x;x],[y;y],[z;z],[col;col],'LineStyle',':' ,'facecol','no', 'edgecol','interp','LineWidth',1  );
    
    surfH{scaf_j}.HitTest= 'off' ;
end
GetHyperB.ScafAllBase = ScafForScatterIndividual; % adapt to multi-scaffold


sH =scatter(ScafForScatterAllSpacing(:,2),5*ScafForScatterAllSpacing(:,1),54 ,'s' ,'filled','Tag','sH') ;
sH.UserData.SFC= ScafForScatterAll ;

% ScafForScatterAll=ScafForScatterOri ;
SFC_C4notation = ScafForScatterAll ;
SFC_C4notation(:,1)= GetHyperB.RelateVec( SFC_C4notation(:,1)) ;
sH.UserData.SFC_C4notation = SFC_C4notation ;
sH.MarkerFaceColor='none' ;
sH.UserData.GetHyperB = GetHyperB ;

%------scaf234556...., closing strand
if strcmp(GetHyperB.UserWantOH,'Yes')
    
    WWW=GetHyperB.ClosingStrand.ExtendScaf ;
    WWW(end+1:end+2,:) =[-1 -1;-1 -1] ;
    %     [~,IndC5] = ismember( WWW(:,1) ,GetHyperB.RelateTable(:,5) ) ;
    %     IndB = GetHyperB.RelateTable(IndC5,1) ;
    %     IndB(GetHyperB.RelateTable(IndC5,2)==-1 )= length(GetHyperB.containBundle)+1 ;
    %     WWW(:,1 ) =  WWW(:,1 )  + (IndB-1)*2 ;
    
    
    CornerNotation=cell(length(GetHyperB.ClosingStrand.ConnScaf2),1) ;   % for closing strand
    Scaf2H=cell(size(CornerNotation) ) ;
    
    
    for k=1:length(GetHyperB.ClosingStrand.ConnScaf2)
        Conn =GetHyperB.ClosingStrand.ConnScaf2{k} ;
        [ ~,ind1] =ismember(Conn(1,:) ,GetHyperB.ClosingStrand.ExtendScaf , 'rows') ;
        [ ~,ind2] =ismember(Conn(2,:) ,GetHyperB.ClosingStrand.ExtendScaf , 'rows')    ;
        if ind1==0 ; ind1= size(WWW,1)-1 ;end
        if ind2==0 ; ind2= size(WWW,1)-1 ;end
        
        if mod(ind1,2)==0 && mod(ind2,2)==0 %not sure
            EE=[WWW(ind2-1,:) ; WWW(ind2,:) ;  WWW(ind1-1,:) ; WWW(ind1,:)];
            
        elseif mod(ind1,2)==1 && mod(ind2,2)==0
            EE=[WWW(ind2-1,:) ; WWW(ind2,:) ;  WWW(ind1,:) ; WWW(ind1+1,:)];
            
        elseif  mod(ind1,2)==0 && mod(ind2,2)==1
            %    EE=[WWW(ind2,:) ; WWW(ind2+1,:) ;  WWW(ind1-1,:) ; WWW(ind1,:)];
            EE=[ WWW(ind1-1,:) ; WWW(ind1,:);  WWW(ind2,:) ; WWW(ind2+1,:) ];   %  ok
        elseif  mod(ind1,2)==1 && mod(ind2,2)==1  %not sure
            EE=[WWW(ind2,:) ; WWW(ind2+1,:) ;  WWW(ind1,:) ; WWW(ind1+1,:)];
            %     EE=[WWW(ind2,:) ; WWW(ind2+1,:) ;  WWW(ind1,:) ; WWW(ind1+1,:)];
            
        end
        EE= EE(EE(:,1)~=-1 ,:) ;
        
        
        %    EE=[WWW(ind2-1,:) ; WWW(ind2,:) ;  WWW(ind1,:) ; WWW(ind1+1,:)];
        CornerNotation{k} = EE  ;
        
        if isempty(EE); continue;end
        
        %-----single extension
        if EE(end,2)==EE(end-1,2)
            EE=EE(1:end-2,:) ;
        elseif EE(1,2)==EE(2,2)
            EE=EE(3:end,:) ;
        end
        %------------
        
        [~,IndC5] = ismember( EE(:,1) ,GetHyperB.RelateTable(:,5) ) ;  % spacing
        IndB = GetHyperB.RelateTable(IndC5,1) ;
        IndB(GetHyperB.RelateTable(IndC5,2)==-1 )= length(GetHyperB.containBundle)+1 ;
        EE(:,1 ) =  EE(:,1 )  + (IndB-1)*2 ;
        %         CornerNotation{k} = EE  ;
        
        if size(unique(EE,'rows') ,1)==1
            EE=[nan nan];
        end
        %         if size(unique(EE,'rows') ,1)~=1
        Scaf2H{k} = plot( EE(:,2),5*EE(:,1) -1,  'k','HitTest','off') ;
        %         else
        
        %         end
        scatter(EE(1,2),5*EE(1,1) -1 ,'ks','filled','HitTest','off');
        
    end
else
    Scaf2H=[];
    CornerNotation=[];
end
GetHyperB.ClosingCornerNotation =CornerNotation ;
axis normal ;
% %----------3D
axes(ax);

surfH3D = cell( length( GetHyperB.scafC5), 1)  ;  % scaffold 3D surf graphic handles in cell
for k = 1 : length(surfH3D)
    surfH3D{k}=GetHyperB.plotScafR_cylindermodelAllBase(1, k) ; surfH3D{k}.HitTest= 'off' ;
end

% surfH3D=GetHyperB.plotScafR_cylindermodel ;

BaseScafR3DXYZ =GetHyperB.plotScafR_BaseByBase(ScafForScatterAll) ;
sH3D=scatter3(BaseScafR3DXYZ(:,1) ,BaseScafR3DXYZ(:,2) ,BaseScafR3DXYZ(:,3),54,'s','filled','Tag','sH3D') ;
sH3D.MarkerFaceColor='none' ;



% [pScaf2,ScafHelix,pScaf_center,ScafBaseCenterHelix ,NVecscaf,scafBundleRout ]=plotScaf2_Helix_V2( GetHyperB,{GetHyperB.scafC5},Isstap ,[0,0,1] ,TM ) ;     % plot scaf strands
Isstap=0 ;  TM=1 ;
[ScafHelixXYZ,ScafBasesCenter  ]=plotScaf2_Helix_V2_noGraphics( GetHyperB,GetHyperB.scafC5,Isstap ,[0,0,1] ,TM ) ;     % get scaf strands coordinate
% return
for k = 1 : length(surfH3D)
    BVec = ScafHelixXYZ{k}-ScafBasesCenter{k} ; d_Bvec= sqrt(BVec(:,1).^2+BVec(:,2).^2+BVec(:,3).^2 ) ; RR = mean(d_Bvec) ;
    BVec = BVec/RR;
    % h_quiver = quiver3( surfH3D.XData(1,:)', surfH3D.YData(1,:)', surfH3D.ZData(1,:)',BVec(:,1),BVec(:,2),BVec(:,3)) ;
    
    surfH3D{k}.UserData.CylRep=[surfH3D{k}.XData(1,:)',surfH3D{k}.YData(1,:)',surfH3D{k}.ZData(1,:)'] ;
    surfH3D{k}.UserData.BVec=BVec ;
end
JsonTxt_scaf.UserData.Mode=1 ;
JsonTxt_scaf.ButtonDownFcn = @(src,evn)JsonTxt_scafButtonDown(src,evn,surfH3D,sH3D) ;
% return

% pStapleH=plotStaple_Helix( GetHyperB  ) ;
Isstap=1 ;  TM=1 ;
[pStapleH,StapHelix,pStap_center,StapBaseCenterHelix ,NVecstap,stapBundleRout ]=plotScaf2_Helix_V2( GetHyperB,GetHyperB.StapList3,Isstap ,[1,0,0] ,TM ) ; % plot 3D staple strands

CornerNotation2=cell(size(CornerNotation)) ;
if strcmp(GetHyperB.UserWantOH,'Yes')
    for k=1:length( CornerNotation)
        EE=    CornerNotation{k} ;
        if isempty(EE);continue;end
        
        if EE(end,2)==EE(end-1,2)
            EE=EE(1:end-2,:) ;
            
        elseif EE(1,2)==EE(2,2)
            EE=EE(3:end,:) ;
        end
        if size(unique(EE,'rows') ,1)==1
            EE=[];
        end
        CornerNotation2{k}=EE ;
    end
    pScaf2H=plotScaf2_Helix( GetHyperB ,CornerNotation2 ) ;
else
    pScaf2H=[];
end

% BaseScafR3DXYZ =GetHyperB.plotScafR_BaseByBase(ScafForScatter) ;
% sH3D=scatter3(BaseScafR3DXYZ(:,1) ,BaseScafR3DXYZ(:,2) ,BaseScafR3DXYZ(:,3),54,'s','filled','Tag','sH3D') ;
% sH3D.MarkerFaceColor='none' ;

for k=1:length(pStapleH)   % 3D staple
    
    %hard code for MEO
    if length( pStapleH{k}.XData)==26 
        plotH{k}.Color= [ 0 0 1 ]  ;
    elseif length( pStapleH{k}.XData)==36
        plotH{k}.Color= [ 1 1 0 ] ;
    end
    
    
    pStapleH{k}.Color =   plotH{k}.Color ;
    delete(pStap_center{k}) ;
end

for k=1:length(pStapleH)   % 3D staple
    pStapleH{k}.UserData.HelicalRep_XYZ = [pStapleH{k}.XData',  pStapleH{k}.YData',pStapleH{k}.ZData' ];
    BVec = StapHelix{k}-StapBaseCenterHelix{k} ; d_Bvec= sqrt(BVec(:,1).^2+BVec(:,2).^2+BVec(:,3).^2 ) ; RR = mean(d_Bvec) ;
    BVec = BVec/RR;
    
    pStapleH{k}.UserData.BVec=BVec ;
    pStapleH{k}.HitTest = 'off' ;
end
JsonTxt_stap.UserData.Mode=1 ;
JsonTxt_stap.ButtonDownFcn = @(src,evn)JsonTxt_stapButtonDown(src,evn,pStapleH) ;

%--------------3D



jsonSlider1.Callback=@(src,evn)scaffold2DChange(src,evn,surfH,surfH3D) ;
jsonSlider2.Callback=@(src,evn)staple2DChange(src,evn,plotH,pStapleH) ;
jsonPop.Callback=@(src,evn)labelChange(src,evn,HeadTail) ;
sH.ButtonDownFcn =@(src,evn)HighLightBase(src,evn,sH3D) ;
sH3D.ButtonDownFcn =@(src,evn)HighLightBase3D(src,evn,sH) ;
ax2.UserData.PlottedScafR = GetHyperB.ScafRouting ;
jsonPop.Value = 4;
labelChange(jsonPop,'',HeadTail) ;
jsonSlider2.Value = 0.4 ;
staple2DChange(jsonSlider2,'',plotH,pStapleH) ;  %-------staple
jsonSlider1.Value=0.9 ;
scaffold2DChange(jsonSlider1,'',surfH,surfH3D) ;

btn2.Callback= @(src,evn)exportjson(src,evn,GetHyperB) ;
btn3.Callback= @(src,evn)exportCSV(src,evn) ;

Str={'--'};
for k=1:length(GetHyperB.scafC5) ;  Str{k}=strcat('Scaffold-',num2str(k),'(', num2str(size(ScafForScatterIndividual{k} ,1) ),')' );    end
Str{k+1}= 'Show all scaffolds' ;
jsonPop2.String = Str; jsonPop2.Value = k+1;
jsonPop2.Callback=@(src,evn)ShowHideScaf(src,evn,surfH3D,surfH) ;




uistack(sH,'top');
uistack(sH3D,'top');


% toc
% fprintf('End of test , \n   ' )
% profile viewer
% return ; %

GetHyperB.pSeq= [] ;
GetHyperB.pSeqExt= [];
for scaf_j = 1 : length( GetHyperB.scafC5)
    
    for k = 1 : length(surfH3D)
        surfH3D{k}.Visible ='off';
        surfH{k}.Visible ='off';
    end
    surfH3D{scaf_j}.Visible ='on';
    surfH{scaf_j}.Visible ='on';
    
    %     UIControl_FontSize_bak = get(0, 'DefaultUIControlFontSize');  %
    %     change default fontsize and turn it back
    %     set(0, 'DefaultUIControlFontSize', 18);
    %     set(0, 'DefaultUIControlFontSize', UIControl_FontSize_bak);
    %
    
    list = {'p7249','p7560','p7704', 'p8064', 'pCS3_L_7560', 'pCS4_7557', 'pCS5_7559' ,'other' ,'randomScaf'};
    strShow = strcat('Select scaffold(',num2str(scaf_j) ,'/', num2str(length( GetHyperB.scafC5)),')',' Current Scaf length =', num2str(size(ScafForScatterIndividual{scaf_j},1 )),'. Containing non-ATCG cause random Scaf seq '  ) ;
    fprintf(strcat(' Current Scaf length =', num2str(size(ScafForScatterIndividual{scaf_j},1 )),'\n') ) ;
    
    opt.Resize = 'on';
    opt.Interpreter='tex' ;
    %     str = strcat('Current long Scaf ~= ' , num2str(NBaseOri) ,' bp' ) ;
    %     str= ['\fontsize{10}' str newline  'Specify the numbers of piceses to cut: ' ];
    
    UIControl_FontSize_bak = get(0, 'DefaultUIControlFontSize');
    set(0, 'DefaultUIControlFontSize', 14);
    % insdati = menu('Can you please help me?','Yes','No')
    
    opts.Interpreter = 'tex';
    [indx,tf] = listdlg('PromptString',strShow,'ListString',list,'SelectionMode','single' ,'ListSize', [500,250] );
    set(0, 'DefaultUIControlFontSize', UIControl_FontSize_bak);
    
    if tf ==1
        switch indx
            case 1
                pSeq = p7249 ;
            case 2
                pSeq = p7560 ;
            case 3
                pSeq = p7704 ;
            case 4
                pSeq = p8064 ;
                
            case 5
                pSeq = pCS3_L_7560 ;
            case 6
                pSeq = pCS4_7557 ;
            case 7
                pSeq = pCS5_7559 ;
                
            case 8
                pSeq = inputdlg('Enter custom sequence:','scaffold', [1 50]) ;
                pSeq= pSeq{1};
                isATCG = or(or(or(ismember(pSeq,'A') , ismember(pSeq,'T')) , ismember(pSeq,'C') ),ismember(pSeq,'G'))  ;
                if sum(isATCG)~=length(isATCG)
                    pSeq= randseq(length(isATCG)) ;
                    fprintf( 'Entering wrong seqence, use random sequence instead  \n')
                end
            case 9
                pSeq=   randseq(size(ScafForScatterIndividual{scaf_j},1)) ;
                
        end
        GetHyperB.pSeq{scaf_j}=pSeq(1 : size(ScafForScatterIndividual{scaf_j},1 ));
    end
    
    if length(pSeq)< size(ScafForScatterIndividual{scaf_j},1 )
        pSeqExt=[pSeq, randseq( size(ScafForScatterIndividual{scaf_j},1 )-length(pSeq) ) ] ;   %
        GetHyperB.pSeqExt{scaf_j}=pSeqExt;
        %     GetHyperB.pSeq=pSeqExt;
    end
    
end  % muti scaffold

for k = 1 : length(surfH3D)
    surfH3D{k}.Visible ='on';
    surfH{k}.Visible ='on';
end

% pSeq;
%     strOO='<html><body bgcolor="#HHHHHH">Hello</body></html>' ;

% SCafBundleR_withSkip = zeros(size(ScafForScatter)) ;
% for k=1:size(SCafBundleR_withSkip,1)
%
%      GetHyperB.RelateTable
% end
% [~,ik] =ismember(ScafForScatter(:,1) ,GetHyperB.RelateTable(:,5)) ;
% SCafBundleR_withSkip=GetHyperB.RelateTable(ik,1) ;



pSeqAll= GetHyperB.pSeqAll ;
ScafForScatter; ScafForScatterAll;
StapSeq = cell( length( GetHyperB.StapList3) ,6 ) ;
StapAllBaaeCell=  cell( length( GetHyperB.StapList3) ,1 ) ;
ScafUnUsed= 1:size(ScafForScatterAll,1)   ;   % ssScaf includes connections and scaffold loops

StapOnWhichScaf =zeros(length( GetHyperB.StapList3) , length( GetHyperB.scafC5) )  ; % matrix form

[MultiScafLength,~  ]=cellfun(@size,ScafForScatterIndividual) ;
accMultiScafLengths = cumsum(MultiScafLength) ;
accMultiScafLength1=[1; accMultiScafLengths(1:end-1) ];


btn4.Enable='off';
for k= 1:length( GetHyperB.StapList3)
    CornerRout = GetHyperB.StapList3{k} ;
    BaseR_notyet =interpolateBase(  CornerRout )  ;
    StapBaseR= setdiff(BaseR_notyet , skipBase,'rows' ,'stable') ;  %C5 notation
    %-----------Introduce insert
    StapBaseRExtra = [StapBaseR , [1: size(StapBaseR,1)]']; 
    insertBaseP = intersect(insertBase,StapBaseR,'rows' ) ;
    [~,bb2] = ismember(insertBaseP,StapBaseR ,'rows'   ) ;
    
    insertBaseExtra = [insertBaseP , bb2+0.5]; 
    Merge= [StapBaseRExtra; insertBaseExtra] ;
    sortMerge = sortrows(Merge,3) ;
    StapBaseR2=sortMerge(:,1:2) ;
    StapBaseR=StapBaseR2;           
   
     %   ---
%     [tf,indbase ] = ismember(StapBaseR2 ,ScafForScatterAll, 'rows' ) ;  %C5
    [tf,indbase ] = ismembertol(StapBaseR2 ,ScafForScatterAll, 'ByRows',true ,'OutputAllIndices',true) ;  %C5
    for sk = 1:length(indbase)-1
        if length(indbase{sk})==1 &&  length(indbase{sk+1})>1
           if indbase{sk}>indbase{sk+1}(1)
            indbase{sk+1} = sort( indbase{sk+1},'descend') ;
           else
            indbase{sk+1} = sort( indbase{sk+1},'ascend') ;   
           end
        end   
    end
    
    [C,ia,ic] = unique(StapBaseR2,'stable','rows') ;
    indbase= cell2mat(indbase(ia) ) ;
    %------
    
%     indbase(ia)
    
    ScafUnUsed=setdiff(ScafUnUsed,indbase) ;
    if sum( tf==0)>0
        fprintf(' has single stranded staple bases, sequence unknown, strand %i  ' ,k)
        btn4.Enable='on';
        indbaseOri =indbase ;
        indbase(indbase==0)=1 ;
        if max(indbase)> length(pSeqAll)
            QQ=seqcomplement(pSeqExt(indbase) ) ;
            QQ(tf==0)='?' ;
            QQ(indbase>length(pSeqAll)) = '*' ;
            fprintf(' scaffold length is not enough, strand %i \n  ' ,k)
        else
            CylindesForSSstap  = StapBaseR(:,1) ;
            [~ ,bCyl] = ismember(CylindesForSSstap ,  GetHyperB.RelateTable(:,5) ) ;
            %             GetHyperB.RelateTable(:,2)~=-1
            IsFakeCyls=   GetHyperB.RelateTable(bCyl,2)==-1;
            QQ=seqcomplement(pSeqAll(indbase) ) ;
            
            QQ(tf==0)='?' ;
            QQ(and(IsFakeCyls==0,tf==0)   ) ='T' ;       % know it's located on cylinders from cylinder model, not pseudo cylinder due to staple overhang.
            
            if sum(QQ=='?')>0
                fprintf(' Detect non-paired staple bases on pseudo-cylinder(Not from cylinder model) . May due to overhang.  ' )
            else
                fprintf(' Detect non-paired staple bases on cylinder model. May due to PolyT (Force them to T).  ' )
            end
        end
        fprintf('\n')
    elseif  max(indbase)>length(pSeqAll)
        
        QQ=seqcomplement(pSeqExt(indbase) ) ;
        QQ(indbase>length(pSeqAll)) = '*' ;
        fprintf(' scaffold length is not enough, strand %i \n  ' ,k)
        
    else
        QQ=seqcomplement(pSeqAll(indbase) ) ;
        
        
        
    end
    %     k
    %     StapSeq{k,3} =  seqcomplement(pSeq(indbase) ) ;
    StapSeq{k,3} = QQ ;
    
    StartP =[GetHyperB.RelateVec(CornerRout(1,1)) ,  CornerRout(1,2) ] ;
    EndP =[GetHyperB.RelateVec(CornerRout(end,1)) ,  CornerRout(end,2) ] ;
    
    StapSeq{k,1} = strcat(num2str(StartP(1)),'[', num2str(StartP(2)),']'  )  ;
    StapSeq{k,2} = strcat(num2str(EndP(1)),'[', num2str(EndP(2)),']'  )  ;
    StapSeq{k,4} = num2str(length(QQ)) ;
    
    %     rgbZero2One =  rand(1,3) ;   % staple color code
    rgbZero2One =  pStapleH{k}.Color ;   % staple color code
    
    Scale256 = rgbZero2One* 255 ;
    HexCode=strcat( dec2hex(round(Scale256(1)),2) , dec2hex(round(Scale256(2)),2) , dec2hex(round(Scale256(3)),2 )   ) ;
    %     StapSeq{k,5} = num2str(size(StapBaseR,1)) ;   % color
    StapSeq{k,5} =HexCode;
    t_json.BackgroundColor(k,:) = rgbZero2One ;
    StapAllBaaeCell{k} =StapBaseR ;
    
    %     BelongScaf = indbase'>accMultiScafLength(1) ;
    
    BelongScaf = sum(and( (indbase-accMultiScafLength1')>0, (indbase-accMultiScafLengths')<0)   )>0 ;
    %     if k== length( GetHyperB.StapList3)-1
    %         fsdf=23
    %     end
    
    StapOnWhichScaf(k,BelongScaf ) =1 ;
    StapSeq{k,6} = strcat('Scaf : ',' ' ,num2str(find( BelongScaf ) )  ) ;
    
end
% % % % % GetHyperB.ScafAllBase = ScafForScatterIndividual; % adapt to
% % % % % multi-scaffold % move to forward for better debugging

GetHyperB.StapAllBase = StapAllBaaeCell;
GetHyperB.ScafUnUsed = ScafUnUsed;
QQ=GetHyperB.ScafUnUsed ;



if ~isempty(QQ) % if have saved overhangs sequences , July 3 2019
    dQQ =diff(QQ)~=1;
    Head=[QQ(1),QQ(find(dQQ)+1 )] ;
    Tails=   [QQ(find(dQQ) ) , QQ(end)];
    HeadTail=[Head;Tails] ;
    UnUsedSeq= cell(size(HeadTail,2) ,1) ;
    for k=1: length(UnUsedSeq)
        UnUsedSeq{k}=pSeqAll(HeadTail(1,k):HeadTail(2,k)) ;
    end
    GetHyperB.ScafUnUsedSeq = UnUsedSeq ;
end

LL= cellfun(@length,StapAllBaaeCell) ;
% %------------ hard code
% StapL_Thres = 10 ;
% if sum(LL<=StapL_Thres) >0
% fprintf('ignore staples with lengths below than %s. Better to run again Preview \n', num2str(StapL_Thres)) ;
% StapAllBaaeCell=StapAllBaaeCell(LL>StapL_Thres) ;
% GetHyperB.StapList3= GetHyperB.StapList3(LL>StapL_Thres) ;
% GetHyperB.HeadOfStap= GetHyperB.HeadOfStap(LL>StapL_Thres,: ) ;
% end
% %---------------


btn4.Callback= @(src,evn)overhangSeq(src,evn,StapSeq,CornerNotation,GetHyperB,ScafForScatterAll,Scaf2H,pScaf2H) ;
% btn3.Cal
t_json.Data =StapSeq;
t_json.UserData.BeforeOHSeq =StapSeq ;
t_json.UserData.TablePrevClick = [-1,-1];
t_json.CellSelectionCallback = @(src,evn)tabelSelectHighLight(src,evn,pStapleH,plotH,pScaf2H,Scaf2H, jsonSlider2,h_bingraph ) ;

% t_json. ButtonDownFcn=@(sr,evn)tabelButton(src,evn) ;
% t_json.CellEditCallback= @(src,evn)tabelButton(src,evn,pStapleH,plotH) ;

fprintf('finish  initial cadnano \n') ;
toc

% setGOfontsize( gctab , 10 , {'UIControl'} )  % set fontsize for uicontrol in this tab
f25=figure(25);clf;
f25.Name='Staple length distribution' ; set(f25,'NumberTitle','off');
histogram(LL,1:1:max(LL));
title(strcat('Staple length distribution. # of strands=' ,num2str(length(LL)), ' Tol. Base=',num2str(sum(LL))  ) ) ;

IndTooLong =find(LL> 60) ;
for k=1:length(IndTooLong)
    fprintf('Staple %i longer than 60\n',IndTooLong(k)) ;
end
if isempty(IndTooLong)
    fprintf('All staple shorter than 60. \n')
end

IndTooShort =find(LL<30 ) ;
for k=1:length(IndTooShort)
    fprintf('Staple %i shorter than 30 (%i).\n',IndTooShort(k) ,LL(IndTooShort(k))) ;
end
if isempty(IndTooShort)
    fprintf('All staple longer than 30 . \n')
end
% uistack(sH,'top');
% uistack(sH3D,'top');

% figure(98);clf; hold on ;Int=50;
% for iB= 1 :max(SCafBundleR_withSkip)
%     subplot(2,max(SCafBundleR_withSkip),iB)
% histogram(find(SCafBundleR_withSkip==iB) ,0:Int:length(SCafBundleR_withSkip)+1 ,'Normalization','probability') ;
% end
% subplot(2,max(SCafBundleR_withSkip),max(SCafBundleR_withSkip)+1:2*max(SCafBundleR_withSkip));hold on
% for iB= 1 :max(SCafBundleR_withSkip)
%
% histogram(find(SCafBundleR_withSkip==iB) ,0:Int:length(SCafBundleR_withSkip)+1 ,'Normalization','probability') ;
% end
% title(sprintf('Scaffold Bases distribution versus bundles , Int= %s ',num2str(Int))); xlabel('Scaffold base index ')

axes(ax);

%------------------------Instruction and legend
hLg= legend([surfH3D{1},pStapleH{1}],'1','2','Location','northwest' ) ;
hLg.String={'\bfscaffold (index \bf\color[rgb]{0.24,0.15,0.66}Start 5'' \color{black}to \color[rgb]{0.78,0.78,0.18}End 3''\rm\color{black})','\bfstaples (other colors) '};
hLg.Interpreter='tex';        %latex
hLg.Orientation='vertical';
%             ForLegend.Marker='.' ; ForLegend.Marker='none';
hLg.ButtonDownFcn=@(src,evn)LegendBoxing_cadnano( src,evn,ax );
hLg.Title.String='Click me for instructions' ;
hLg.Units='normalized'; %hLg.AutoUpdate ='off';
hLg.Position=[0.0063 0.9528 0.1569 0.0387];
%------------------------
if UseStapleGraph ==1
    h_bingraph.ButtonDownFcn=@(src,evn)StapGraphButtonDown(src,evn,t_json,  pStapleH,  plotH ,jsonSlider2) ;
end
% profile viewer
% opts.Interpreter = 'tex';
% opts.Default = 'No';
% answer = questdlg('\fontsize{15} Inspect Scaffold and Staple routing mapping ?'  , ...
%     'Inspect or not?', ...
%     'Yes','No',opts);
% if strcmp(answer ,'Yes')
%     GetHyperB.InspectRouting(t_json,  pStapleH,  plotH ,jsonSlider2);
% end

fprintf('End of preview  \n') ;

end

function ShowHideScaf(src,evn,surfH3D,surfH)
for k = 1 : length(surfH3D)
    surfH3D{k}.Visible ='off';
    surfH{k}.Visible ='off';
end
if src.Value >length(surfH3D)
    for k = 1 : length(surfH3D)
        surfH3D{k}.Visible ='on';
        surfH{k}.Visible ='on';
    end
else
    surfH3D{ src.Value}.Visible ='on';
    surfH{ src.Value}.Visible ='on';
end
end


function overhangSeq(src,evn,StapSeq,CornerNotation,GetHyperB,ScafForScatter,Scaf2H,pScaf2H)
t_json=findobj(gcf,'Tag','t_json') ;
t_OH=findobj(gcf,'Tag','OHTable');  CompareCol=[2,5,6,7,8] ; IsSame=true;
FirstRow = {t_OH.Data{1,2}, t_OH.Data{1,5},  t_OH.Data{1,6}, t_OH.Data{1,7}, t_OH.Data{1,8}} ;
for coli=2:size(t_OH.Data,1)
    NextRow =  {t_OH.Data{coli,2}, t_OH.Data{coli,5},  t_OH.Data{coli,6}, t_OH.Data{coli,7}, t_OH.Data{coli,8}} ;
    for k=1:5
        if ismember(k,[1,2,5])
            TF=strcmp(FirstRow{k} ,NextRow{k} ) ;  %str
        else
            TF=FirstRow{k}==NextRow{k};     % num
        end
        
        IsSame=  and(IsSame,TF) ;
        if ~TF; break;end
    end
end


opts.Interpreter = 'tex';
opts.Default = 'random generate';
skipBase= GetHyperB.skipBase ;

if ~isempty(GetHyperB.SeqsOfClosing)
    str='Closing Strand Seq '  ;
    answerClosingSeq = questdlg('Found Closing Strand Seq. Do you want to use it ?',str, ...
        'Yes','OverWrite' ,'Yes');
else
    answerClosingSeq='OverWrite' ;
end


if strcmp(answerClosingSeq,'OverWrite')
    switch IsSame
        case false
            answer = questdlg({'\fontsize{15} Overhangs detected!! How would you like for overhangs ?','Overhang data is not identical.' }  , ...
                'Overhang options', ...
                'Optimized Seq','Custom assign',opts);
        case true
            answer = questdlg('\fontsize{15} Overhangs detected!! How would you like for overhangs ?', ...
                'Overhang options', ...
                'Optimized Seq','Custom assign','Apply the same sequence',opts);
    end
    
    switch answer
        case 'Optimized Seq'
            SeqsOfClosing = cell(length(CornerNotation),1) ;   % need to be created for all cases
            ssScaf= GetHyperB.ScafUnUsedSeq ;
            EntireScaf=GetHyperB.pSeq ;
            OtherOverhangs=[];
            for k =1: length(CornerNotation)
                Corn=  CornerNotation{k};
                if isempty(Corn) ; continue; end
                
                %-------old, random sequences
                %             BaseByBase =interpolateBase(  Corn )  ;
                %             RandSeq = randseq(size(BaseByBase,1)) ;
                %-----------------------------------
                if  abs(diff(Corn(1:2,2)))+1==1  % single entesion
                    OneOverHang1='';
                else
                    [OneOverHang1,~]=getOverSeq_ssScaf( abs(diff(Corn(1:2,2)))+1 ,OtherOverhangs,ssScaf,EntireScaf) ;
                end
                OtherOverhangs=[string(OtherOverhangs); string(OneOverHang1 )] ;
                if size(Corn,1)==4
                    if abs(diff(Corn(3:4,2)))+1==1  % single entesion
                        OneOverHang2= '';
                    else
                        [OneOverHang2,~]=getOverSeq_ssScaf( abs(diff(Corn(3:4,2)))+1 ,OtherOverhangs,ssScaf,EntireScaf) ;
                    end
                else
                    OneOverHang2='';
                end
                
                OtherOverhangs=[string(OtherOverhangs); string(OneOverHang2 )] ;
                SeqsOfClosing{k} =  [OneOverHang1,OneOverHang2]  ;
                %              SeqsOfClosing{k} = RandSeq   ;
            end
        case 'Apply the same sequence'
            SeqsOfClosing = cell(length(CornerNotation),1) ;   % need to be created for all cases
            Corn=  CornerNotation{1};
            %         BaseByBase =interpolateBase(  Corn )  ;
            %         RandSeq = randseq(size(BaseByBase,1)) ;
            
            ssScaf= GetHyperB.ScafUnUsedSeq ;
            EntireScaf=GetHyperB.pSeq ;
            OtherOverhangs=[];
            [OneOverHang1,~]=getOverSeq_ssScaf( abs(diff(Corn(1:2,2)))+1 ,OtherOverhangs,ssScaf,EntireScaf) ;
            OtherOverhangs=[string(OtherOverhangs); string(OneOverHang1 )] ;
            [OneOverHang2,~]=getOverSeq_ssScaf( abs(diff(Corn(3:4,2)))+1 ,OtherOverhangs,ssScaf,EntireScaf) ;
            OtherOverhangs=[string(OtherOverhangs); string(OneOverHang2 )] ;
            SeqsOfClosing{k} =  [OneOverHang1,OneOverHang2]  ;
            
            for k =1: length(CornerNotation)
                SeqsOfClosing{k}= [OneOverHang1,OneOverHang2];
                %              SeqsOfClosing{k}=RandSeq;
            end
            
        case 'Custom assign' 
            SeqsOfClosing = cell(length(CornerNotation),1) ;   % need to be created for all cases
            opts.Interpreter = 'tex'; opts.Resize='on';
            
            GetHyperB.ClosingStrand.ssT_shift = zeros( length(CornerNotation),1 ) ;
            for k =1: length(CornerNotation)
                Corn=  CornerNotation{k};
                if ~isempty(Corn)
                    if Corn(end,2)==Corn(end-1,2)   % single extension
                        Corn=Corn(1:end-2,:) ;
                    elseif Corn(1,2)==Corn(2,2)
                        Corn=Corn(3:end,:) ;
                    end
                    if size(unique(Corn,'rows') ,1)==1
                        continue;
                    end
                    
                    BaseByBase =interpolateBase(  Corn )  ;
                else
                    BaseByBase=[];
                end
                Scaf2H{k}.Color=[1,0,0]; pScaf2H{k}.Color=[1,0,0];
                Scaf2H{k}.LineWidth =2 ; pScaf2H{k}.LineWidth =2 ;
                str=strcat('\fontsize{12}#OH= ', num2str(k) ,'/',num2str(length(CornerNotation))  ,', Enter Seq,  N=', num2str(size(BaseByBase,1)),...
                    ', highlighted as red, incorrect input will be randomly assigned. Key in ''self-comp or 0-4 (length of ssT )'' for staple-staple overhang complementary.')   ;
                answer = inputdlg(str, 'Sequence', [2 100] , {randseq(size(BaseByBase,1))} ,opts) ;
                ssT_lengths={'0', '1' , '2', '3','4'} ; 
                if ~isempty(answer)
                    if strcmp( answer{1}, 'self-comp' ) || ismember(answer{1}, ssT_lengths )
                        if strcmp( answer{1}, 'self-comp' )
                            ssTL = 0 ;
                        else
                            ssTL = str2num(answer{1}) ;
                        end
                        hL = size(BaseByBase,1)/2;
                        RootCS = GetHyperB.ClosingStrand.RootAndExtC5(:, [1 3]) ;
                        IsHeadAsRoot = ismember(Corn(1,:) , RootCS  ,'rows' ) ;
                        IsTailAsRoot = ismember(Corn(end,:) , RootCS  ,'rows' )  ;
                        
                        if xor(IsHeadAsRoot , IsTailAsRoot)  % assume as double overhang case1
                            fprintf('Assigning staple overhangs in this pair complementary to each other (Double overhang1).\n')
                            halfSeq = randseq(size(BaseByBase,1)/2) ;
                            answer{1}  = [halfSeq , seqrcomplement(halfSeq)] ;
                            GetHyperB.ClosingStrand.ssT_shift(k) = ssTL ;
                            
                            answer{1}(1:ssTL)='A' ;  % force these staple bases to T (CS='A') ;
                            answer{1}(end-ssTL+1:end)='A' ;
                            if ~IsHeadAsRoot
                                answer{1}(1:hL)= circshift( answer{1}(1:hL) , -ssTL) ;
                            else
                                answer{1}(hL+1:end)= circshift( answer{1}(hL+1:end) , ssTL) ;
                            end
                        else   % assume as double overhang case2 
                            fprintf('Assigning staple overhangs in this pair complementary to each other (Double overhang2).\n')
%                             sdfsf=2 ;
                            halfSeq = randseq(size(BaseByBase,1)/2) ;
                            answer{1}  = [halfSeq , seqrcomplement(halfSeq)] ;
                            GetHyperB.ClosingStrand.ssT_shift(k) = ssTL ;
                            if IsHeadAsRoot %Insert 'T' in both ends which are roots 
                            answer{1}(1:ssTL)='A' ;  % force these staple bases to T (CS='A') ;
                            answer{1}(end-ssTL+1:end)='A' ;
                            else  %Insert 'T' in the middle,  which are roots 
                            answer{1}(hL-ssTL+1:hL)='A' ;  % force these staple bases to T (CS='A') ;
                            answer{1}(hL+1:hL+ssTL)='A' ;
                            end

                            
                        end
                        
                    end
                    
                    
                    tf1 =length(answer{1})==size(BaseByBase,1) ;
                    indA=answer{1}=='A' ;  indT=answer{1}=='T' ;  indC=answer{1}=='C' ; indG=answer{1}=='G' ;
                    isATCG=or(or(or(indA, indT) ,indC ), indG) ;
                    tf2=  sum(isATCG==1)==length(isATCG) ;
                    
                    if ~tf1 || ~tf2
                        answer{1}  = randseq(size(BaseByBase,1)) ;
                    end
                else
                    answer{1}  = randseq(size(BaseByBase,1)) ;
                end
                SeqsOfClosing{k} = answer{1}   ;
                Scaf2H{k}.Color=[0,0,0]; pScaf2H{k}.Color=[0,0,0];
                Scaf2H{k}.LineWidth=1 ; pScaf2H{k}.LineWidth=1 ;
            end
    end
    
else
    SeqsOfClosing= GetHyperB.SeqsOfClosing;
end




% StapSeq = cell( length( GetHyperB.StapList3) ,6 ) ;
ScafForScatter ;
ScafForScatterExt=ScafForScatter ;
[aScaf,bScaf] = cellfun(@size, GetHyperB.pSeq) ;
if  length(GetHyperB.pSeqAll)>=size(ScafForScatter,1)
    ExtScafseq = GetHyperB.pSeqAll(1:size(ScafForScatter,1)) ; %removed unused scaf, otherwise indexes mess up.
else  %imported scaffold not enough
    ExtScafseq = [GetHyperB.pSeqAll ,  randseq(size(ScafForScatter,1)-length(GetHyperB.pSeqAll) ) ] ; % add some ramdom scaffold, latter will be subs as *
end

for k=1:length(CornerNotation)   %put closing strand sequence after used scaffold
    EE=    CornerNotation{k} ;
    if ~isempty(EE)
        if EE(end,2)==EE(end-1,2)   % single extension
            EE=EE(1:end-2,:) ;
        elseif EE(1,2)==EE(2,2)
            EE=EE(3:end,:) ;
        end
        if size(unique(EE,'rows') ,1)==1
            EE=[] ;
        end
        BaseByBase =interpolateBase(  EE )  ;
        
    else % incase particularly add nicks by overhangs
        BaseByBase=[];
    end
    ScafForScatterExt=[ScafForScatterExt;BaseByBase];
    ExtScafseq=[ExtScafseq,SeqsOfClosing{k}] ;
end

% ScafForScatter;
StapSeq = t_json.UserData.BeforeOHSeq ;
for k= 1:length(StapSeq)
    k;
    CornerRout = GetHyperB.StapList3{k}  ;
    BaseR_notyet =interpolateBase(  CornerRout )  ;
    StapBaseR= setdiff(BaseR_notyet , skipBase,'rows' ,'stable') ;  %C5 notation
    
    [tf,indbase ] = ismember(StapBaseR ,ScafForScatterExt, 'rows' ) ;
    
    
    if sum( tf==0)>0
        indOri =indbase ;
        indbase(indbase==0) = 1 ;
        QQ=seqcomplement(ExtScafseq(indbase) ) ;
        indExceedScaf = indbase > length(GetHyperB.pSeqAll) ;
        indNotOverHang =  indbase <= size(ScafForScatter,1) ;
        
        QQ(and(indExceedScaf,indNotOverHang )) ='*';
        QQ(indOri==0) ='T';
        %         fprintf('Should not happen !!! \n ')
        fprintf('Found single-stranded staple without scaffold, Not due to overhang. Assign These bases sequece as T (PolyT case). Strand = %i !!! \n ',k) ;
    else
        %         if k==44
        %             sdf=3
        %         end
        
        QQ=seqcomplement(ExtScafseq(indbase) ) ;
        indExceedScaf = indbase > length(GetHyperB.pSeqAll) ;
        indNotOverHang =  indbase <= size(ScafForScatter,1) ;
        QQ(and(indExceedScaf,indNotOverHang )) ='*';
        
    end
    StapSeq{k,3} = QQ ;
end

nCell = size(StapSeq,1) ; t_json.BackgroundColor=t_json.BackgroundColor(1:nCell,:) ;
StapSeq{nCell+1,3}='Below are closing strands ';




t_json.BackgroundColor(nCell+1,:)  =[1,1,1] ;

for k=1 : length(CornerNotation)
    Corn=  CornerNotation{k};
    if isempty(Corn)
        t_json.BackgroundColor(nCell+k+1,:)  =[0.7,0.7,0.7] ;
        continue;
    end
    
    BaseByBase =interpolateBase(  Corn )  ;
    StartP =[GetHyperB.RelateVec(BaseByBase(1,1)) ,  BaseByBase(1,2) ] ;
    EndP =[GetHyperB.RelateVec(BaseByBase(end,1)) ,  BaseByBase(end,2) ] ;
    StapSeq{nCell+1+k,1} = strcat(num2str(StartP(1)),'[', num2str(StartP(2)),']'  )  ;
    StapSeq{nCell+1+k,2} = strcat(num2str(EndP(1)),'[', num2str(EndP(2)),']'  )  ;
    StapSeq{nCell+1+k,3} = SeqsOfClosing{k}  ;
    StapSeq{nCell+1+k,4} = num2str(length(SeqsOfClosing{k}  ));
    t_json.BackgroundColor(nCell+k+1,:)  =[0.7,0.7,0.7] ;
end

t_json.Data =StapSeq;
GetHyperB.pSeqExt=ExtScafseq;
GetHyperB.SeqsOfClosing=SeqsOfClosing ;
fprintf('Finish Overhang Seq. \n')

end

% function mydialog(CornerNotation)
% d = dialog('Unit','normalized','Position',[0.05 0.1 0.2 0.4],'Name','My Dialog');
% movegui(d,'center')
% txt = uicontrol('Parent',d,'Unit','normalized','Style','text',....
%     'Position',[0.1 0.8 0.3 0.1], 'String','Click the close button when you''re done.');
%
% btn1 = uicontrol('Parent',d, 'Unit','normalized','Position',[0.1 0.2 0.1 0.1],....
%     'String','1',  'Callback','delete(gcf)');
% t_Seq = uitable('Parent',d, 'Unit','normalized','Position',[0.3 0.1 0.65 0.8],...
%     'ColumnEditable', true, 'ColumnName',{'Required L' , 'Seqence'});
% for k =1: length(CornerNotation)
%     Corn=  CornerNotation{k};     BaseByBase =interpolateBase(  Corn )  ;
%     t_Seq.Data{k,1} = size(BaseByBase,1)   ;
% end
%
%
% end


function  exportCSV(src,evn)
t_json=findobj(gcf,'Tag','t_json') ;

prompt = {'Enter File name:'};  dlg_title = 'Input';
num_lines = 1;  defaultans = {'CSVfilename'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

fileID = fopen([pwd filesep answer{1} '.csv' ],'w');
fprintf(fileID , 'Staple#,Start,End,Seq,Length,Color,Note \n'  )    ;
for k=1:size(t_json.Data,1)
    if contains(t_json.Data{k,3} , '?')
    fprintf(fileID , '%i,%s,%s,%s,%s,#%s,%s \n' ,k ,t_json.Data{k,1},t_json.Data{k,2},t_json.Data{k,3},t_json.Data{k,4},t_json.Data{k,5}, 'This strand contains unknow Seqs.'  )    ;
    else
    fprintf(fileID , '%i,%s,%s,%s,%s,#%s,%s \n' ,k ,t_json.Data{k,1},t_json.Data{k,2},t_json.Data{k,3},t_json.Data{k,4},t_json.Data{k,5},t_json.Data{k,6})    ;
    end
end
fclose(fileID);

end



% function tabelButton(src,evn)
% %  evn.Indices
% % % src.UserData.prev
% % % selectStrand =  evn.Indices(1)  ;
% %     c = uisetcolor  ;
% %     pStapleH{evn.Indices(1)}.Color = c ;
% %     plotH{evn.Indices(1)}.Color = c ;
% %     src.BackgroundColor(evn.Indices(1),:) =c ;
% % src.Data{evn.Indices(1),evn.Indices(2)}='';
% sdfsdf=3
% end


function   tabelSelectHighLight(src,evn,pStapleH,plotH,pScaf2H,Scaf2H,jsonSlider2 ,h_bingraph )
% src.UserData
% evn.Indices
src.UserData.SelectedCell= evn.Indices ;
src.UserData.pStapleH= pStapleH ;
src.UserData.plotH= plotH ;

if ~isempty(evn.Indices)
    selectStrand =  evn.Indices(1)  ;
    for k=1:length( pStapleH)
        if k== selectStrand %%|| k==126 || k==12
            pStapleH{k}.LineWidth =5;        plotH{k}.LineWidth =5 ;
            pStapleH{k}.Color(4) = 1 ;      plotH{k}.Color(4) = 1 ;
        else
            pStapleH{k}.LineWidth =2 ;       plotH{k}.LineWidth =2 ;
            pStapleH{k}.Color(4) =jsonSlider2.Value ;  plotH{k}.Color(4)=jsonSlider2.Value ;
        end
    end
    
    for k2=1:length( pScaf2H)  % closing strand
        if  isempty( pScaf2H{k2})
            continue ;
        end
        
        if k2+k+1 == selectStrand
            pScaf2H{k2}.LineWidth =4 ;        Scaf2H{k2}.LineWidth =4 ;
        else
            pScaf2H{k2}.LineWidth =1 ;        Scaf2H{k2}.LineWidth =1 ;
        end
        
    end
    src.UserData.prev=selectStrand;
    
    %     if evn.Indices(2)==5  &&   src.UserData.TablePrevClick(1)~=evn.Indices(1)  % select color column
    %         Sti =evn.Indices(1) ;
    %         c = uisetcolor ;
    %         pStapleH{Sti}.Color =c ;
    %         plotH{Sti}.Color =c ;
    %         src.BackgroundColor(Sti,:)=c ;
    %         Scale256 = c* 255 ;
    %         HexCode=strcat( dec2hex(round(Scale256(1)),2) , dec2hex(round(Scale256(2)),2) , dec2hex(round(Scale256(3)),2 )   ) ;
    %
    %         src.Data{evn.Indices(1),5} = HexCode ;
    %         src.UserData.TablePrevClick = evn.Indices ;
    %     end
end
% evn.Indices
if ~isempty(h_bingraph)  && ~isempty(evn.Indices)
    h_bingraph.MarkerSize=3;
    highlight(h_bingraph,evn.Indices(1));
    
    %-----
    % sdfsf=3
    h_bingraph.NodeLabel=[] ; % label cadnano starting point
    edgebinsWeak = conncomp(h_bingraph.UserData.Graph,'OutputForm','cell','Type','weak') ;
    edgesWeak = conncomp(h_bingraph.UserData.Graph,'OutputForm','vector','Type','weak') ;
    Neighbors =  edgebinsWeak{edgesWeak(evn.Indices(1))}  ;
    
    QQ = distances(h_bingraph.UserData.Graph,Neighbors',Neighbors')  ;
    QQ(QQ==Inf)=0 ;
    max(max(QQ)) ;
    [Start,End]=find(QQ==max(max(QQ))) ;
    %     P = shortestpath(h_bingraph.UserData.Graph,Neighbors(Start),Neighbors(End)) ;
    P = shortestpath(h_bingraph.UserData.Graph,Neighbors(Start(1)),Neighbors(End(1))) ;
    
    Neighbors=P ;  % update order
    % labelnode(h_bingraph,1:729,  cellstr(num2str([1:729]')) )
    strs = src.Data(Neighbors,1) ;
    LL = src.Data(Neighbors,4) ;
    for k =1:length(strs)
        strs{k} =strcat( strs{k}, ' L=',  LL{k}) ;
    end
    fprintf('\n')
    for k=1:length(LL)-1
        fprintf('%s -',LL{k}) ;
    end
    fprintf('%s',LL{length(LL)});
    fprintf('\n')
    
    
    for k=1:length(strs)
        fprintf('|%s |',strs{k}) ;
    end
    fprintf('\n')
    
    
    % labelnode(h_bingraph,Neighbors,strs) ;
    
end
% if  evn.Indices(2) ==6
%
% %     sdfs=3
%     c = uisetcolor  ;
%     pStapleH{evn.Indices(1)}.Color = c ;
%     plotH{evn.Indices(1)}.Color = c ;
%     src.BackgroundColor(evn.Indices(1),:) =c ;
%
% end




% sdfs=3 ;
end




% function HighLightBase(src,evn,sH3D)
% evn.IntersectionPoint ;
% XY=[src.XData; src.YData]' ;
% ax= findobj(0,'Tag','json3D') ;
% ax2=findobj(0,'Tag','json2D') ;
%
%
%
% [~,ind] = ismember( evn.IntersectionPoint(1:2) ,XY ,'rows' )  ;
%
%
% if ~isfield( src.UserData , 'Highlight')
%     axes(ax2) ;
%      src.UserData.Highlight= scatter(evn.IntersectionPoint(1) ,evn.IntersectionPoint(2),120,'or' ,'filled') ;
%     axes(ax) ;
%      src.UserData.Highlight3D= scatter3(sH3D.XData(ind) ,sH3D.YData(ind),sH3D.ZData(ind),120,'or' ,'filled') ;
%      src.UserData.ind=ind ;
%
%
% else
%      src.UserData.Highlight.XData =evn.IntersectionPoint(1) ;
%      src.UserData.Highlight.YData =evn.IntersectionPoint(2) ;
%      src.UserData.Highlight3D.XData = sH3D.XData(ind) ;
%      src.UserData.Highlight3D.YData = sH3D.YData(ind) ;
%      src.UserData.Highlight3D.ZData = sH3D.ZData(ind) ;
%      src.UserData.ind=ind ;
% end
%
%
% end


function exportjson(src,evn,GetHyperB)
GetHyperB.ExportJSON;
end

function labelChange(src,evn,HeadTail)
k=1;
switch src.Value
    case 1   % show num
        set(HeadTail{1}, 'Visible','on') ;  set(HeadTail{4}, 'Visible','on') ;
        HeadTail{2}.Visible ='off' ;
        HeadTail{3}.Visible ='off' ;
    case 2 % show ends
        set(HeadTail{1}, 'Visible','off') ;  set(HeadTail{4}, 'Visible','off') ;
        HeadTail{2}.Visible ='on' ;
        HeadTail{3}.Visible ='on' ;
    case 3 % both
        set(HeadTail{1}, 'Visible','on') ;  set(HeadTail{4}, 'Visible','on') ;
        HeadTail{2}.Visible ='on' ;
        HeadTail{3}.Visible ='on' ;
    case 4 % none
        set(HeadTail{1}, 'Visible','off') ;  set(HeadTail{4}, 'Visible','off') ;
        HeadTail{2}.Visible ='off' ;
        HeadTail{3}.Visible ='off' ;
end
end


function scaffold2DChange(src,evn,surfH,surfH3D)
for k = 1 : length(surfH)
    surfH{k}.EdgeAlpha=src.Value ;
    surfH3D{k}.EdgeAlpha=src.Value ;
end
end

function staple2DChange(src,evn,plotH,pStapleH)
for k=1: length(plotH)
    plotH{k}.Color(4) = src.Value ;
    pStapleH{k}.Color(4) = src.Value ;
end
end

function JsonTxt_scafButtonDown(src,~,surfH3D,sH3D)
% evn
% src
XYZAll = [];
if src.UserData.Mode == 1   % cyliner rep to helical
    for k=1: length(surfH3D)
        XYZ= surfH3D{k}.UserData.CylRep +  0.8*surfH3D{k}.UserData.BVec ;
        surfH3D{k}.XData=[XYZ(:,1)' ;XYZ(:,1)' ];
        surfH3D{k}.YData=[XYZ(:,2)' ;XYZ(:,2)' ];
        surfH3D{k}.ZData=[XYZ(:,3)' ;XYZ(:,3)' ];
        src.UserData.Mode=2 ;
        XYZAll=[XYZAll;XYZ];
    end
else                            % helical rep to cylinder
    for k=1: length(surfH3D)
        XYZ= surfH3D{k}.UserData.CylRep  ;
        surfH3D{k}.XData=[XYZ(:,1)' ;XYZ(:,1)' ];
        surfH3D{k}.YData=[XYZ(:,2)' ;XYZ(:,2)' ];
        surfH3D{k}.ZData=[XYZ(:,3)' ;XYZ(:,3)' ];
        src.UserData.Mode=1 ;
        XYZAll=[XYZAll;XYZ];
    end
end
sH3D.XData= XYZAll(:,1)'  ;
sH3D.YData= XYZAll(:,2)'  ;
sH3D.ZData= XYZAll(:,3)'  ;


end

function JsonTxt_stapButtonDown(src,evn,pStapleH)
if src.UserData.Mode == 1   % helical rep to cylinder
    for k=1: length(pStapleH)
        XYZ = pStapleH{k}.UserData.HelicalRep_XYZ -  pStapleH{k}.UserData.BVec ;
        pStapleH{k}.XData= XYZ(:,1)' ;
        pStapleH{k}.YData= XYZ(:,2)' ;
        pStapleH{k}.ZData= XYZ(:,3)' ;
    end
    src.UserData.Mode=2 ;
    
else   % cyliner rep to helical
    for k=1: length(pStapleH)
        XYZ = pStapleH{k}.UserData.HelicalRep_XYZ -  0.2*pStapleH{k}.UserData.BVec  ;
        pStapleH{k}.XData= XYZ(:,1)' ;
        pStapleH{k}.YData= XYZ(:,2)' ;
        pStapleH{k}.ZData= XYZ(:,3)' ;
    end
    
    src.UserData.Mode=1 ;
end
end

function StapGraphButtonDown(src,evn,table,  pStapleH,  plotH ,jsonSlider2)

d= [src.XData' ,src.YData' ] -  ones(length(src.YData),1)*evn.IntersectionPoint(1:2) ;
d= d(:,1).^2+d(:,2).^2 ;
Ind= find(d==min(d));

edgebinsWeak = conncomp(src.UserData.Graph,'OutputForm','cell','Type','weak') ;
edgesWeak = conncomp(src.UserData.Graph,'OutputForm','vector','Type','weak') ;

Neighbors =  edgebinsWeak{edgesWeak(Ind)}  ;

EdgeList = table2array(src.UserData.Graph.Edges) ;

IndEdge1 = ismember( EdgeList(:,1),Neighbors) ; %find(IndEdge1)
% IndEdge2 = ismember( EdgeList(:,2),Neighbors)
% evn
if evn.Button==1
    highlight(src,EdgeList((IndEdge1),1),EdgeList((IndEdge1),2) ,'LineWidth',5  )  ;
elseif evn.Button==3
    highlight(src,EdgeList((IndEdge1),1),EdgeList((IndEdge1),2) ,'LineWidth',2  )  ;
    
    return
end



Neighbors =  edgebinsWeak{edgesWeak(Ind)}  ;

for k=1:length(pStapleH)
    if ismember(k,Neighbors )
        pStapleH{k}.LineWidth =5;        plotH{k}.LineWidth =5 ;
        pStapleH{k}.Color(4) = 1 ;      plotH{k}.Color(4) = 1 ;
    else
        pStapleH{k}.LineWidth =2 ;       plotH{k}.LineWidth =2 ;
        pStapleH{k}.Color(4) =jsonSlider2.Value ;  plotH{k}.Color(4)=jsonSlider2.Value ;
    end
end
QQ = distances(src.UserData.Graph,Neighbors',Neighbors') ;
QQ(QQ==Inf)=0 ;
max(max(QQ)) ;
[Start,End]=find(QQ==max(max(QQ))) ;
if length(Start)>1
    Start=Start(1) ;    End=End(1) ;
end

P = shortestpath(src.UserData.Graph,Neighbors(Start),Neighbors(End)) ;
global stapleLength
Neighbors =P ;
strs = table.Data(Neighbors,1) ;
LL = table.Data(Neighbors,4) ;
for k =1:length(strs)
    strs{k} =strcat( strs{k}, ' L=',  LL{k}) ;
end
fprintf('\n');  sum( cellfun(@str2num ,LL)) ;
% fprintf('\n'); sum(cell2mat(LL));
fprintf('Total length =%i in %i staples\n',sum(cellfun(@str2num,LL)), length(LL) ) ;
stapleLength =cellfun(@str2num ,LL) ;

for k=1:length(LL)-1
    fprintf('%s -',LL{k}) ;
end
fprintf('%s',LL{length(LL)}) ;
fprintf('\n')

% for k=1:length(LL)-1
% fprintf('(%s)%s -',num2str(k),LL{k}) ;
% end
% fprintf('(%s)%s',num2str(length(LL)),LL{length(LL)}) ;
% fprintf('\n')

for k=1:length(strs)
    fprintf('|%s |',strs{k}) ;
end
fprintf('\n')

getGlobalstapleLength;
end




