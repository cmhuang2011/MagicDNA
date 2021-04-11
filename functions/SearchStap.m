function  SearchStap(src,evn )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% profile on  % last optimized, 08/26/2019

fH=gcf ;
fprintf('\n') ;
ax= findobj(gcf,'Tag','MechStaple3D');
axes(ax); cltab ; axes(ax); drawnow;
ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;
% GetHyperB.findRT;

% tt=findobj(gcf,'Tag','OHTable') ;
tic
fprintf(' start looking for staple routing Again\n')
fbar = waitbar(0,'looking for staple routing...'); pause(.5) ;

ss_Staple= findobj(fH,'Tag','ss_Staple');
StapOption= ss_Staple.UserData.StapOption ;

%---------
if StapOption.choice.type.straightuncut.Value==1
    StapOption.type='straightuncut' ;
elseif StapOption.choice.type.zigzaguncut.Value==1
    StapOption.type='zigzaguncut' ;
elseif StapOption.choice.type.zigzagcut.Value==1
    StapOption.type='zigzagcut' ;
elseif StapOption.choice.type.straightcut.Value==1
    StapOption.type='straightcut' ;
end

if StapOption.choice.halfXover.yes.Value==1
    StapOption.halfXover='yes';
else
    StapOption.halfXover='no';
end
%---------------------------
GetHyperB.StapOption= StapOption ;  %update to object

GetHyperB.RelateTable=GetHyperB.RelateTableOrr  ;

[~,index]=sortrows(GetHyperB.RelateTable,5);
Vec=GetHyperB.RelateTable(:,4);
TTRelateVec=Vec(index);
WQER =union(TTRelateVec,[],'stable');
GetHyperB.RelateVec=WQER';
GetHyperB.CustomSkipAndInsertion = false; % whenever using routing algorithm, reset to Default
for k =1:length(GetHyperB.containBundle)
GetHyperB.containBundle{k}.skipPosition=[];
GetHyperB.containBundle{k}.InsertPosition=[];    
end

%-----------
% GetHyperB.StapList3=GetHyperB.StapList2;

% return


%
%
waitbar(.2,fbar,'getting initial path....'); figure(fH);

type=2;
GetHyperB.FindStap(type );    % Get properties: StapList ......., stapBP ,stapBPinCell
GetHyperB.ConnectStapleIf_0ntOnScaf ;  % optional, if ssDNA scaffold assigned as 0, input&output: StapList




switch StapOption.type
    case 'straightuncut'
        GetHyperB.StapList3 =GetHyperB.StapList;
    case 'zigzaguncut'
        waitbar(.3,fbar,'Applying staple Xovers....');
        GetHyperB.FindStapStep2;           %Get property:StapList2
        GetHyperB.StapList3 = GetHyperB.StapList2 ;
    case 'zigzagcut'
        waitbar(.3,fbar,'Applying staple Xovers....');
        GetHyperB.FindStapStep2;           %Get property:StapList2
        waitbar(.5,fbar,'Cutting staple 1/2....');
        
        %          GetHyperB.ApartStaples;     %Get property:StapList3
        GetHyperB.ApartStaplesSimple;     %Get property:StapList3
        GetHyperB.StapList2= GetHyperB.StapList3;
        waitbar(.6,fbar,'Cutting staple 2/2....');
        GetHyperB.ApartStaplesSimple([],[]) ;     %Get property:StapList3
        
        
        
    case 'straightcut'
        GetHyperB.StapList2 = GetHyperB.StapList ;
        
        
        % GetHyperB.StapList3= GetHyperB.StapList2;
        GetHyperB.ApartStaplesSimple;     %Get property:StapList3
        GetHyperB.StapList2= GetHyperB.StapList3;
        waitbar(.5,fbar,'Cutting staple....');
        GetHyperB.ApartStaplesSimple([],[]) ;     %Get property:StapList3
        GetHyperB.StapList3 = CalibStapDir( GetHyperB.StapList3,GetHyperB.RelateVec);
        
        %------
        GetHyperB.ConnectStaple3If_0ntOnScaf ;  %  if ssDNA scaffold assigned as 0,input&output: StapList3
        GetHyperB.ApartStaplesSimple([],[]) ;     %Get property:StapList3
        GetHyperB.StapList3 = CalibStapDir( GetHyperB.StapList3,GetHyperB.RelateVec);
        
end

% stepAssign=3 ;Ex_DecoratePlatesWithOverhangs ;

% GetHyperB.ConvertStap('Square');    % Get properties: DigitStapSQ,   HeadOfStep
% GetHyperB.ConvertStap('Honeycomb');    % Get properties: DigitStapHC
% 
% 
% plotH=DrawStapp3(GetHyperB.StapList3,GetHyperB,1,findobj(fH,'Tag','MechStaple2D'),[]) ;
% return


%             output:StapleCell%   Hard Code

% GetHyperB.FindStapStep2;           %Get property:StapList2



% GetHyperB.ApartStaples;     %Get property:StapList3
% GetHyperB.StapList3=GetHyperB.StapList2;
% GetHyperB.StapList3=GetHyperB.StapList;



% figure ; [plotH,HeadTail]=DrawStapp3(GetHyperB.StapList3,[],1,gca,[]) ;
% set(HeadTail{1},'Visible','off') ;
% set(HeadTail{4},'Visible','off') ;
if strcmp(GetHyperB.UserWantOH,'Yes')
    waitbar(.7,fbar,'Extending overhangs....');
    GetHyperB.extendOverhang(ax);   %Rtable chnages
       
end


% figure ; [plotH,HeadTail]=DrawStapp3(GetHyperB.StapList3,[],1,gca,[]) ;
% set(HeadTail{1},'Visible','off') ;
% set(HeadTail{4},'Visible','off') ;

%-------overwrite closing extension, 
% ConvertScafG(GetHyperB);
% ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
% ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC


GetHyperB.ConvertStap('Square');    % Get properties: DigitStapSQ,   HeadOfStep
GetHyperB.ConvertStap('Honeycomb');    % Get properties: DigitStapHC

ConvertScafG(GetHyperB);
ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC

fprintf(' found staple routing \n')
waitbar(.9,fbar,'found staple routing ....');

% toc


% GetHyperB.ConvertStap('Square');    % Get properties: DigitStapSQ,   HeadOfStep
% GetHyperB.ConvertStap('Honeycomb');    % Get properties: DigitStapHC

% GetHyperB.ExportJSON;

% figure(3123);clf;  xx=gca;
% DrawStapp3(GetHyperB.StapList3,[],1,xx,[]) ;

plotH=DrawStapp3(GetHyperB.StapList3,GetHyperB,1,findobj(fH,'Tag','MechStaple2D'),[]) ;
waitbar(.95,fbar,'Visualizing results...');




% return ;

%-----------start plot scaffold
axes(ax);

%
hold on; axis  equal;
Isstap = 0; TM=1 ;
[pScaf2,ScafHelix,pScaf_center,ScafBaseCenterHelix ,NVecscaf,scafBundleRout ]=plotScaf2_Helix_V2( GetHyperB,GetHyperB.scafC5,Isstap ,[0,0,1], TM  ) ;     % plot scaf strands
for k= 1: length(pScaf_center)
    pScaf_center{k}.Visible='off' ;
end
%-----------------------start to plot staples
Isstap=1 ;
[pstap3,StapHelix,pStap_center,StapBaseCenterHelix ,NVecstap,stapBundleRout ]=plotScaf2_Helix_V2( GetHyperB,GetHyperB.StapList3,Isstap ,[1,0,0], TM  ) ;     % plot scaf strands
for k=1:length(pstap3)
    pstap3{k}.Color=  plotH{k}.Color ;
    delete(pStap_center{k}) ;
end
% plotStaple_Helix( GetHyperB  );

toc
%     profile viewer
close(fbar) ;

hLg= legend([pScaf2{1},pstap3{1}],'x','Click me for instructions','Location','northwest' ) ; hLg.String={'scaffold(s) ','staples(other colors) '};
hLg.Interpreter='tex';        %latex
hLg.Orientation='horizontal';
%             ForLegend.Marker='.' ; ForLegend.Marker='none';
hLg.ButtonDownFcn=@(src,evn)LegendBoxing_Stap( src,evn,ax );
hLg.Title.String='Click me for instructions' ;
hLg.Units='normalized'; %hLg.AutoUpdate ='off';
hLg.Position=[0.0063 0.9528 0.1569 0.0387];




% profile viewer



end

