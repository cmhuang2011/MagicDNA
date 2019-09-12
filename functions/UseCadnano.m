function UseCadnano(src,evn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

sdfsf=2;

ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;
ss_Staple=findobj(gcf,'Tag','ss_Staple') ;
ax_scaf= findobj(gcf,'Tag','MechScaffold3D');
axes(ax_scaf); cltab ; axes(ax_scaf);

JsonEx = jsonObject; 
fprintf('Loading cadnano Json file!! \n')  ;


ss_Staple.UserData.JsonEx=JsonEx ;


% if strcmp(GetHyperB.UserWantOH,'Yes')
%    GetHyperB.RelateTable=   GetHyperB.RelateTableOrr ;
% end

% JsonEx=ss_Staple.UserData.JsonEx ;


StapCornerC4= JsonEx.StapCornerRep ;   % still C4 Rep
ScafCornerC4= JsonEx.ScafCornerRep ;   % still C4 Rep

%------convert stap C4 to C5
NewStpList= cell(size(StapCornerC4)) ;
for k=1:length(StapCornerC4)
    CornerC4 =StapCornerC4{k} ;
    CornerC5=CornerC4 ;
    [a,b]=ismember(CornerC4(:,1),GetHyperB.RelateTable(:,4) ) ;    
    CornerC5(:,1) = GetHyperB.RelateTable(b,5) ;
    NewStpList{k} =CornerC5;
end
%--------------
%------convert scaf C4 to C12

CornerC4 =ScafCornerC4{1} ;
CornerC5= [zeros(size(CornerC4,1),1) ,CornerC4 ];
[a,b]=ismember(CornerC4(:,1),GetHyperB.RelateTable(:,4) ) ;
CornerC5(:,1:2) = GetHyperB.RelateTable(b,1:2) ;
NewScaf =CornerC5;
GetHyperB.ScafRouting =  NewScaf ;   % update
GetHyperB.Scaf_fromJSON= NewScaf ;
ConvertScafG(GetHyperB);
ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC
GetHyperB.plotScafR_cylindermodel(1) ;       
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
    GetHyperB.LoadJsonRecoverClosingStrand;   
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

