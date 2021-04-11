function UseCadnano(src,evn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

sdfsf=2;
fprintf('start of UseCadnano\n')  ;

ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;
ss_Staple=findobj(gcf,'Tag','ss_Staple') ;
ax_scaf= findobj(gcf,'Tag','MechScaffold3D');
axes(ax_scaf); cltab ; axes(ax_scaf);

JsonEx = jsonObject;
fprintf('Loading cadnano Json file!! \n')  ;

answer = questdlg('Would you update skips and insertions in caDNAno?', ...
    'skip/insertion', ...
    'Yes','No','Yes');
if strcmp(answer,'Yes')
    GetHyperB.CustomSkipAndInsertion = true ;
    AllskipC4 = JsonEx.skipC4 ;   
%     AllskipC4=zeros(0,2); %hard
    AllinsertC4 = JsonEx.insertC4 ; 
    
    %-----skip
    AllskipBCB_C5=zeros(size(AllskipC4,1),4) ;  % [Bundle Cylinder C5 base]
    [a,b]=ismember(AllskipC4(:,1),GetHyperB.RelateTable(:,4) ) ;
    AllskipBCB_C5(:,1:2) = GetHyperB.RelateTable(b,1:2) ;
    
    AllskipBCB_C5(:,3) = GetHyperB.RelateTable(b,5) ;
    AllskipBCB_C5(:,4) = AllskipC4(:,2) ;
    for k=1: max(AllskipBCB_C5(:,1))
        GetHyperB.containBundle{k}.skipPosition =  AllskipBCB_C5(AllskipBCB_C5(:,1)==k , :) ;
    end
    %-----insert
    AllinsertBCB_C5=zeros(size(AllinsertC4,1),4) ;  % [Bundle Cylinder C5 base]
    [a,b]=ismember(AllinsertC4(:,1),GetHyperB.RelateTable(:,4) ) ;
    AllinsertBCB_C5(:,1:2) = GetHyperB.RelateTable(b,1:2) ;
    
    AllinsertBCB_C5(:,3) = GetHyperB.RelateTable(b,5) ;
    AllinsertBCB_C5(:,4) = AllinsertC4(:,2) ;
    for k=1: max(AllinsertBCB_C5(:,1))
        GetHyperB.containBundle{k}.InsertPosition =  AllinsertBCB_C5(AllinsertBCB_C5(:,1)==k , :) ;
    end
    
else
    GetHyperB.CustomSkipAndInsertion = false ;    % use default
end

ss_Staple.UserData.JsonEx=JsonEx ;


% if strcmp(GetHyperB.UserWantOH,'Yes')
%    GetHyperB.RelateTable=   GetHyperB.RelateTableOrr ;
% end

% JsonEx=ss_Staple.UserData.JsonEx ;


StapCornerC4= JsonEx.StapCornerRep ;   % still C4 Rep
ScafCornerC4= JsonEx.ScafCornerRep ;   % still C4 Rep

% load('ScafCOrnerMOE_mod.mat','ScafCOrnerMOE_mod') ; % hard code
% load('StapCOrnerMOE_mod.mat','StapCOrnerMOE_mod') ; % hard code
% % load('ScafCOrnerMOE.mat','ScafCOrnerMOE') ; % hard code
% % load('StapCOrnerMOE.mat','StapCOrnerMOE') ; % hard code
% ScafCornerC4 = ScafCOrnerMOE_mod ;
% StapCornerC4 = StapCOrnerMOE_mod ;

%------convert stap C4 to C5
NewStpList= cell(size(StapCornerC4)) ;
for k=1:length(StapCornerC4)
    CornerC4 =StapCornerC4{k} ;
    CornerC5=CornerC4 ;
    [a,b]=ismember(CornerC4(:,1),GetHyperB.RelateTable(:,4) ) ;
%         b
%         if sum(b==0)>1
%             sdfs=32
%         end
        
    CornerC5(:,1) = GetHyperB.RelateTable(b,5) ;
    % Errors are possible if users call algorithm again after update caDNAno routing
    % with overhangs, Due to RelateTable has been reset without overhang
    % cylinder. In that case, recall staple routing algorithm to get
    % RelateTable with overhang cylinder.
    
    NewStpList{k} =CornerC5;
end
%--------------
%------convert scaf C4 to C12

% [aa,~]=cellfun(@size,ScafCornerC4 ) ;

% if sum(aa<= 6)>0
%     ExtendScaf=[];
%     UpdateClosingStrand = ScafCornerC4(aa<= 6)  ;
%     for k= 1:length(UpdateClosingStrand)
%     ExtendScaf=[ExtendScaf;UpdateClosingStrand{k}];
%     end
%     GetHyperB.ClosingStrand.ExtendScaf = ExtendScaf ; % update for cadnanoinitial,
% end

% dsaf=3
% ScafCornerC4=ScafCornerC4(aa> 6)    ;   % hard code to exclude modification of closing strand

GetHyperB.ScafRouting= cell(length(ScafCornerC4),1) ;
GetHyperB.Scaf_fromJSON= cell(length(ScafCornerC4),1) ;

for k = 1 : length(ScafCornerC4)
    CornerC4 =ScafCornerC4{k} ;
    CornerC5= [zeros(size(CornerC4,1),1) ,CornerC4 ];
    [a,b]=ismember(CornerC4(:,1),GetHyperB.RelateTable(:,4) ) ;
    CornerC5(:,1:2) = GetHyperB.RelateTable(b,1:2) ;
    NewScaf =CornerC5;
    GetHyperB.ScafRouting{k} =  NewScaf ;   % update
    GetHyperB.Scaf_fromJSON{k}= NewScaf ;
end

ConvertScafG(GetHyperB);
ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC
% GetHyperB.plotScafR_cylindermodel(1) ;
% NewScaf=
%---------------

% NewStpList =
GetHyperB.SaveStapList3 =GetHyperB.StapList3 ;

GetHyperB.StapList3=NewStpList ;
% sdf=3
GetHyperB.ConvertStap('Square');    % Get properties: DigitStapSQ,   HeadOfStep
GetHyperB.ConvertStap('Honeycomb');    % Get properties: DigitStapHC


if strcmp(GetHyperB.UserWantOH,'Yes')
    %    GetHyperB.RelateTable=   GetHyperB.RelateTableOrr ;
%     GetHyperB.LoadJsonRecoverClosingStrand;
end



JsonEx.checksandwich ;

%-------------hard code to check staple lengths distr.
% Bundles=[7, 8 ,9 ; 10,12 ,11 ] ;
% StapleLengths = cell(size(Bundles)) ;
% for k=1:length(NewStpList)
%      CornerRout =NewStpList{k} ;
%
%     Bi =  GetHyperB.RelateTable(CornerRout(1,1)==GetHyperB.RelateTable(:,5) ,1 ) ;
%
%     BaseR_notyet =interpolateBase(  CornerRout )  ;
%     if  ismember(Bi ,Bundles )
%         [u,v]= find(Bi==Bundles) ;
%         StapleLengths{u,v} = union( StapleLengths{u,v}  , size(BaseR_notyet,1)) ;
%
%     end
% end
% figure ; clf ;
% for k=1:2
%     for k2 =1 :3
%    subplot(2,3, 3*(k-1)+k2 ) ;
%    histogram( StapleLengths{k,k2}  );
%    title(strcat('Bundle ',num2str(Bundles( k,k2)) , '{ }', 'Avg=',num2str(mean(StapleLengths{k,k2}))) ) ;
%    fprintf('Bundle %i L: %s \n',Bundles( k,k2), num2str(StapleLengths{k,k2})     )
%     end
% end
%--------------




fprintf('Scaffold routing is also updated!! \n')  ;

fprintf('end of UseCadnano\n')  ;
end

