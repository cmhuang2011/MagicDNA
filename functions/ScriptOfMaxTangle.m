% script of maximize the minimum the number of connecting staples between
% different scaffold. For multi-purpose use, heuristic opt.
tic
ss_Assembly= findobj(0,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;
 GetHyperB.Scaf_fromJSON = GetHyperB.Scaf_fromCadDOM  ;

ScafXovers=GetHyperB.getXoverinScaf( GetHyperB.Scaf_fromCadDOM) ;
CheckBasePosition = ScafXovers(:,3:3:12) ;
Ns = zeros(size(CheckBasePosition,1) ,1 ) ;
for k=1: length(Ns)
 Ns(k) = length( unique(CheckBasePosition(k,:))) ;   
 unique(CheckBasePosition(k,:));
 if abs(diff (  unique(CheckBasePosition(k,:))))~=1
   Ns(k) = 5 ;
 end
end
ScafXovers=ScafXovers(Ns==2,:) ;   % avoid force Xovers
% return
% C4ScafOXover = zeros(size(ScafXovers,1) , 8) ;
% C4ScafOXover(:,2:2:8 ) =ScafXovers(:,3:3:12);
% for k= 1:size(C4ScafOXover,1)
%     A1=  ScafXovers(k,1:2) ;
%     [AA,BB] = ismember(A1, GetHyperB.RelateTable(:,1:2) ,'rows') ;
%     C4ScafOXover(k,1)= GetHyperB.RelateTable(BB,4) ;
%     
%     A2=  ScafXovers(k,4:5) ;
%     [AA,BB] = ismember(A2, GetHyperB.RelateTable(:,1:2) ,'rows') ;
%     C4ScafOXover(k,3)= GetHyperB.RelateTable(BB,4) ;
%     
%     A3=  ScafXovers(k,7:8) ;
%     [AA,BB] = ismember(A3, GetHyperB.RelateTable(:,1:2) ,'rows') ;
%     C4ScafOXover(k,5)= GetHyperB.RelateTable(BB,4) ;
%     
%     A4=  ScafXovers(k,10:11) ;
%     [AA,BB] = ismember(A4, GetHyperB.RelateTable(:,1:2) ,'rows') ;
%     C4ScafOXover(k,7)= GetHyperB.RelateTable(BB,4) ;
% end

NOriGPairBCB= GetHyperB.Scaf_fromJSON ;


% TwoCycles=    removeScafXover(GetHyperB,NOriGPairBCB{1},Xover)  ;
for k = 1:size(ScafXovers,1)
Xover= [ScafXovers(k,1: 6) ;ScafXovers(k,7:12) ];
NOriGPairBCB=    removeScafXover_general(GetHyperB,NOriGPairBCB,Xover)  ;
end
%                 NOriGPairBCB= GetHyperB.AddBridgeXover2TwoCycle(NOriGPairBCB,Xover,ScafXovers) ;

GetHyperB.ScafRouting =NOriGPairBCB ;
% return
%----------
% Output=IntegrateScaffold(GetHyperB ,NOriGPairBCB) ;
% [ Goodcc,whichCylcc ] = CheckScafRXover(Output{1}, 6,[] ) 

% figure(22);clf; ax=gca;
%  GetHyperB.drawBCBdata(Output,ax);
% GetHyperB.ScafRouting =Output;
% ConvertScafG(GetHyperB);
% ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
% ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC

% Output=IntegrateScaffold(GetHyperB ,GetHyperB.ScafRouting) ;
% [ Goodcc,whichCylcc ] = CheckScafRXover(Output{1}, 8,[] ) 

%----------
% return
N =1 ;
Results = zeros(N, 6) ;
RoutingCell = cell(N,1) ;
for cc= 1:N
    
% NOriGPairBCB = GetHyperB.ScafRouting;
% if length(NOriGPairBCB) ==1
%     [OneXover, XoverList]=Given2CylceFindXoverList(GetHyperB,NOriGPairBCB{1},NOriGPairBCB{1},[]) ;
%     XoverList=XoverList( randperm(size(XoverList,1) ) ,:) ;
%     for k = 1:1 %:  size(XoverList,1)
%         AllRouting =NOriGPairBCB ;
%         OneXoverFromList= [XoverList(k,1:6) ; XoverList(k,7:12)] ;
%         AllRoutingNew=GetHyperB.AddXover2OneCycle(AllRouting ,OneXoverFromList) ;
%         
%         OneXover=GetHyperB.Given2CylceFindXoverList(AllRoutingNew{1},AllRoutingNew{2},[]) ;
%         
%         if isempty(OneXover)
%         break;
%             
%         end
%         [NOriGPairBCB,UnUsedX, CycleInvolved] = GetHyperB.AddXover22cycles(AllRoutingNew,OneXover)
%         
%         
%     end
%     
% end

Output=IntegrateScaffold(GetHyperB ,NOriGPairBCB) ;
GetHyperB.Scaf_fromCadDOM =Output;
fprintf(' Got a new long scaffold !!\n')
return

[OutputScafTangle,CutGood] = SplitScafToMultiScaf()   ;  %% input = GetHyperB.Scaf_fromCadDOM ;

if CutGood==1
Results(cc,:) = [ cc , OutputScafTangle ,min(OutputScafTangle)] ;
RoutingCell{cc} =GetHyperB.ScafRouting ;
cc
end

end
% SearchScaf

toc
save MaxMinConnect1000.mat Results RoutingCell
