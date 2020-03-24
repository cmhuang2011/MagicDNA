% script of maximize the minimum the number of connecting staples between
% different scaffold. For multi-purpose use, heuristic opt.
tic;clc ;
fH=gcf;
%%  Find out all internal Xover which MagicDNA found out before to split the entire scaffold routings to 70~80 cycles, remaining the manual-modified features at Tail.
% Later the integration of so many cycles is a random process which
% enhence the heuristic optimization.

ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ; % get the handle of the class and access to all data.

ScafXovers=GetHyperB.getXoverinScaf( GetHyperB.Scaf_fromJSON) ;  % start with the routing which has been modified through caDNAno.
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
Scaf_Cycles= GetHyperB.Scaf_fromJSON ;

for k = 1:size(ScafXovers,1)
    Xover= [ScafXovers(k,1: 6) ;ScafXovers(k,7:12) ];
    Scaf_Cycles=    removeScafXover_general(GetHyperB,Scaf_Cycles,Xover)  ;
end

GetHyperB.ScafRouting =Scaf_Cycles ;

%----visualize the currting scaffold routing if needed.
fShow = figure ; clf ; subplot(1,3,1) ;GetHyperB.plotScafR_cylindermodelMulti(1 ,'IsoColor') ;  % 1: current   ,2: from MagicDNA/CadDOM

% return
%% This second portion of this script involves {the integration, split, evaluation} in a loop to find out the 4 scaffold routings with maximized tangle.
fprintf('Doing second part---------------------- \n')
%----------
% return
figure(fH) ;
N =10 ;
Results = zeros(N, 6) ;
RoutingCell = cell(N,1) ;
    

for cc= 1:N
    
    
    Output_OneCycle=IntegrateScaffold(GetHyperB ,Scaf_Cycles) ;    %  the integration
    GetHyperB.Scaf_fromCadDOM =Output_OneCycle;
    %-
%     figure(fShow); subplot(1,3,2) ; GetHyperB.plotScafR_cylindermodelMulti(2 ,'IsoColor') ; figure(fH) ;
    fprintf(' Got a new long scaffold !!\n')
    % return
    %--------------------
    [OutputScafTangle,CutGood] = SplitScafToMultiScaf_noQuery()   ;  %% split in to 4, evaluation
    
    if CutGood==1
        Results(cc,:) = [ cc , OutputScafTangle ,min(OutputScafTangle)] ;
        RoutingCell{cc} =GetHyperB.ScafRouting ;
        cc
        
        %-
%         figure(fShow); subplot(1,3,3) ; GetHyperB.plotScafR_cylindermodelMulti(1 ,'IsoColor') ;  figure(fH) ;
    end
    
end
% SearchScaf

toc
% save MaxMinConnect1000.mat Results RoutingCell

%%  Visualize the objective value.
return 
figure ; plot(Results(:,6) )