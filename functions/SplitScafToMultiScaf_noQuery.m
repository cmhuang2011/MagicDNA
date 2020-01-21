function [OutputScafTangle, CutGood] =SplitScafToMultiScaf_noQuery(src,evn)
%UNTITLED Summary of this function goes here
%   multi-scaffold 2nd function: spliting one single loop into multiple with length constraint. 
% profile on
% tic
ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;

%     prompt = { ['\fontsize{10}Enter the one Scaffold to split:' newline ' 1->MagicDNA       2->caDNAno'] };   
%     dlgtitle = 'Input';
%     dims = [1 50];
%     definput = {'1'};
%     opt.Resize = 'on';
%     opt.Interpreter='tex' ;
%     answer = inputdlg(prompt,dlgtitle,dims,definput,opt)  ; 
% if strcmp(answer{1},'1')
% GetHyperB.ScafRouting =  GetHyperB.Scaf_fromCadDOM ;  %  use routing from MagicDNA. 
% % MagicDNA=1
% else
% GetHyperB.ScafRouting = GetHyperB.Scaf_fromJSON ;  %  use routing from caDNAno. 
% % cadnano=1
% end

GetHyperB.ScafRouting =  GetHyperB.Scaf_fromCadDOM ;  %  use routing from MagicDNA. 


if length(GetHyperB.ScafRouting) ==1   % start with one cycle case
    fprintf('\n')
%     fprintf('doing something .....\n')
    fprintf('-----------------------------------------------------\n')
    NOriGPairBCB = GetHyperB.ScafRouting;
    BaseRoutOri = interpolateBase_ThreeCol( NOriGPairBCB{1} ); NBaseOri = size(BaseRoutOri,1) ;
    str = strcat('Current long Scaf ~= ' , num2str(NBaseOri) ,' bp' ) ;
    str= ['\fontsize{10}' str newline  'Specify the numbers of piceses to cut: ' ];
    
%     prompt = { str,'\fontsize{10}Enter accuracy (%):'};   
%     dlgtitle = 'Input';
%     dims = [ 1 50 ; 1 50 ];
%     definput = {'3','10'};
%     opt.Resize = 'on';
%     opt.Interpreter='tex' ;
%     answer = inputdlg(prompt,dlgtitle,dims,definput,opt)  ; Npiece = str2num(answer{1}) ; Acc = str2num(answer{2}) ;    
    %----if want to call this function without query N and % for heauristic optmization, commned the previous line and assign values to the two variables. 
        Npiece=4 ; Acc=5 ;

%     return
    [~, XoverList]=Given2CylceFindXoverList(GetHyperB,NOriGPairBCB{1},NOriGPairBCB{1},[]) ;
    
    IndAllOri =QuerryXoverInd_BCB( XoverList,BaseRoutOri ) ;
    Closeness = abs( IndAllOri(:,1)-IndAllOri(:,2) ) ;
    
    %-----------

    %     figure; hist(Closeness)
    Threshold = NBaseOri/Npiece*0.8 ;
    XoverIndRemain=  and(Closeness> Threshold,Closeness< (NBaseOri-Threshold)) ;
    XoverList= XoverList(XoverIndRemain ,:) ;
%     XoverList= XoverList(Closeness> Threshold ,:) ;
%     XoverList= XoverList(Closeness< (NBaseOri-Threshold) ,:) ;
    fprintf('Exclude Xovers which bridge two positions on the one scaffold within %i bases. \n',Threshold)
    
    
    %-------------for animation
%      IndAllOri= IndAllOri(Closeness> Threshold ,:) ;   
%     skipBase= GetHyperB.skipBase ;
%     [~,MapCyl] = ismember(skipBase(:,1) ,GetHyperB.RelateTable(:,5)) ;
%     skipBaseBCB =[GetHyperB.RelateTable(MapCyl,1:2) ,skipBase(:,2) ] ;
%     BaseRoutOriWSkip= setdiff(BaseRoutOri , skipBaseBCB,'rows' ,'stable') ;  %BCB notation
%     IndAllOri =QuerryXoverInd_BCB( XoverList,BaseRoutOriWSkip ) ;
%     save('IndAllOri.mat','IndAllOri')
    %------------
    
    
    
    XoverList=XoverList( randperm(size(XoverList,1) ) ,:) ;  % (random sequence).
    fprintf('check applying %i Xovers to split the scaf routing and match the range ( %i percent) .....\n', size(XoverList,1),Acc )
    fprintf('NBaseOri~= %i  Mean: %i  Lower: %i .\n' ,NBaseOri ,round((NBaseOri/Npiece)) ,round((NBaseOri/Npiece)*(1-Acc/100)) ) ;
%     fprintf('NBaseOri~= %i  Mean: %i  Lower: %i Higher: %i .\n' ,NBaseOri ,round((NBaseOri/Npiece)) ,round((NBaseOri/Npiece)*(1-Acc/100)) ,round((NBaseOri/Npiece)*(1+Acc/100) )) ;
    
    AllRouting =NOriGPairBCB ;
    allGood= zeros( Npiece-1 ,1) ;
    for Nc = 1 : Npiece-1  % outer loop for cutting times
        for k = 1:1 :  size(XoverList,1)  % inner loop for applying all possible Xovers(random sequence).
            OneXoverFromList= [XoverList(k,1:6) ; XoverList(k,7:12)] ;
            AllRoutingNew=GetHyperB.AddXover2OneCycle(AllRouting ,OneXoverFromList);
            if length(AllRoutingNew) == Nc+1  % make sure the number of total cylce is continously increasing.
                
                CheckAll =zeros( length(AllRoutingNew) ,1 ) ; % lengths of each scaffold
                for cc= 1 : length(AllRoutingNew)
                    BaseRout2 = interpolateBase_ThreeCol( AllRoutingNew{cc} ); CheckAll(cc) = size(BaseRout2,1) ;
                end
                CheckAll  ;% lengths after splitting
                
                %----------here to define the length constraint. 
                Good =   and( CheckAll>(NBaseOri/Npiece)*(1-Acc/100) , CheckAll<(NBaseOri/Npiece) ) ;
%                 Good =   and( CheckAll>(NBaseOri/Npiece)*(1-Acc/100) ,  CheckAll<(NBaseOri/Npiece)*(1+Acc/100) ) ;
                
                % example of customized legnth rule : 
%               %  Good= [CheckAll(1)<7557,CheckAll(2)<8064] ; % hard code
                if sum(Good)>=Nc
                    OneXoverFromList ;
                    IndS =QuerryXoverInd_BCB( OneXoverFromList,BaseRoutOri ) ;
                    %                     fprintf('Ind %i : %s \n' ,Nc ,num2str(IndS(:)'));                    
                    fprintf('Cut %i : %s \n' ,Nc ,num2str(CheckAll'));
                    allGood(Nc) = 1;
                    AllRouting= AllRoutingNew;

                    break ;
                end
            end
        end
    end
    CutGood =0 ;
    if sum(allGood)== length(allGood)  % update if every cutting is good
        CutGood=1;
        GetHyperB.ScafRouting=AllRouting ;
        ConvertScafG(GetHyperB);
        ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
        ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC
        fprintf('Found the crossover satisfying the ranges, update data .....\n') ;
    else
        %             figure ; ax=gca;
        %              GetHyperB.drawBCBdata2(AllRouting,157)
        %              drawBCBdata2(obj,cycleList,num)
        OutputScafTangle=[];
        fprintf('Cutting not good \n')
        return;
    end
    
    if ~isempty( GetHyperB.StapAllBase)  % evaluate cropping scaffold with staples
%         fprintf('Checking sacffold and staple mapping. \n')
        OutputScafTangle= GetHyperB.InspectRouting_export;
%      sdfsf=3
    end
else
    
    fprintf('The number of loops is not one. Nothing happens. \n')
end

% profile viewer
fprintf('Finish spliting Scaffold ...... \n')
% toc
end

