function SplitScafToMultiScaf2(src,evn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Give up, Use V0




ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;  

 GetHyperB.ScafRouting =  GetHyperB.Scaf_fromCadDOM ;
if length(GetHyperB.ScafRouting) ==1   % start with one cycle case
    fprintf('\n')
    fprintf('doing something .....\n')
    
    NOriGPairBCB = GetHyperB.ScafRouting;
    BaseRoutOri = interpolateBase_ThreeCol( NOriGPairBCB{1} ); NBaseOri = size(BaseRoutOri,1) ;
    
    
    str = strcat('Current long Scaf ~= ' , num2str(NBaseOri) ,' bp' ) ;
    str= ['\fontsize{10}' str newline  'Specify the numbers of piceses to cut: ' ];
    
    prompt = {  str,'\fontsize{10}Enter accuracy (%):'};
    dlgtitle = 'Input';
    dims = [1 50 ; 1 50];
    definput = {'4','10'};
    opt.Resize = 'on';
    opt.Interpreter='tex' ;
    answer = inputdlg(prompt,dlgtitle,dims,definput,opt)  
    Npiece = str2num(answer{1}) ; Acc = str2num(answer{2}) ;
    
    %     return
     fprintf('NBaseOri~= %i  Mean: %i  Lower: %i Higher: %i .\n' ,NBaseOri ,round((NBaseOri/Npiece)) ,round((NBaseOri/Npiece)*(1-Acc/100)) ,round((NBaseOri/Npiece)*(1+Acc/100) )) ;
   
    
    [OneXover, XoverList]=Given2CylceFindXoverList(GetHyperB,NOriGPairBCB{1},NOriGPairBCB{1},[]) ;
    
    
%      XoverList=XoverList( randperm(size(XoverList,1) ) ,:) ;  
%      IndsXoverCorners= zeros(size(BaseRoutOri,1)  ,4 ) ;
     [~,b1] =ismember(XoverList(:,1:3),  BaseRoutOri,'rows' );
     [~,b2] =ismember(XoverList(:,4:6),  BaseRoutOri,'rows' );
     [~,b3] =ismember(XoverList(:,7:9),  BaseRoutOri,'rows' );
     [~,b4] =ismember(XoverList(:,10:12),  BaseRoutOri,'rows' );
     IndsXoverCorners= [b1,b2,b3,b4]  ;
     XoverCrossRange = abs( IndsXoverCorners(:,4) -IndsXoverCorners(:,2) )  ;
     
     FilterXover1 = and( XoverCrossRange>(NBaseOri/Npiece)*(1-Acc/100) , XoverCrossRange<(NBaseOri/Npiece)*(1+Acc/100))  ; 
     FilterXover2 = and( (NBaseOri-XoverCrossRange)>(NBaseOri/Npiece)*(1-Acc/100) ,  (NBaseOri-XoverCrossRange)<(NBaseOri/Npiece)*(1+Acc/100))  ; 
    
     FilterXover=or(FilterXover1,FilterXover2) ; sum(FilterXover)
     XoverList2 = XoverList(FilterXover ,:)  ;
%      XoverCrossRange2 = XoverCrossRange(FilterXover ,:) ;
     IndsXoverCorners2= IndsXoverCorners(FilterXover ,:)  ;
     
     IndRedundant = IndsXoverCorners2(:,1)< IndsXoverCorners2(:,2) ;
     XoverList2= XoverList2(IndRedundant ,:) ;
     IndsXoverCorners2= IndsXoverCorners2(IndRedundant ,:) ;
   
%      IndsXoverCorners2(IndSwitch, :) =  IndsXoverCorners2(IndSwitch, [3 4  1 2 ]) ; 
%      XoverList2(IndSwitch , :) = XoverList2(IndSwitch , [7:9 10:12 1:3 4:6   ] ) ;
     
     
     [a,b] = sortrows(IndsXoverCorners2) ;
     IndsXoverCorners2=IndsXoverCorners2(b,:)  ;
     XoverList2=XoverList2(b,:)  ;
        
     IndCrossingStartPoint_Use =diff(IndsXoverCorners2(:,1:2),1,2) ;
     IndCrossingStartPoint=   ~and( IndCrossingStartPoint_Use>(NBaseOri/Npiece)*(1-Acc/100) , IndCrossingStartPoint_Use<(NBaseOri/Npiece)*(1+Acc/100))  ; 
        XoverList2= [XoverList2(IndCrossingStartPoint==0,:) ; XoverList2(IndCrossingStartPoint==1,:)  ]  ;
        
        CrossOnesNeedFlipInd =  IndsXoverCorners2(IndCrossingStartPoint==1,:)  ;
        CrossOnesNeedFlipInd=CrossOnesNeedFlipInd(:,[3,4,1,2]) ; CrossOnesNeedFlipInd(:,2)= CrossOnesNeedFlipInd(:,2)+NBaseOri;
        IndsXoverCorners2= [IndsXoverCorners2(IndCrossingStartPoint==0,:) ; CrossOnesNeedFlipInd  ];
 
        
                [ 1:size(IndsXoverCorners2,1) ; IndsXoverCorners2']' 
                [ 1:size(XoverList2,1) ; XoverList2']' 
                

     sdsf=3 ;
     tic
     ChoseXover= nchoosek(1:size(IndsXoverCorners2,1),Npiece-1 ) ;
     JudgeChose = zeros(size(ChoseXover,1) ,1 ) ;
     for k =1:length(JudgeChose)
        IndsOnScaf = IndsXoverCorners2(  ChoseXover(k, :) ,1:2)  ;
        DiffIndsOnScaf= diff([IndsOnScaf(2:end ,1) ,IndsOnScaf(1:end-1,2) ]  ,1 ,2)  ;       
        if sum(DiffIndsOnScaf<0) ==length(DiffIndsOnScaf)
        JudgeChose(k)= 1 ;
        end
     end
     toc
      fprintf(' number of  No Crossing=  %i  %i \n' , sum(IndCrossingStartPoint==0) ,sum(IndCrossingStartPoint==1) );
    
%      size(QQ)
     QQ =ChoseXover(JudgeChose==1,:) ;   
     fprintf(' number of  allowable solution =  %i \n' ,size(QQ,1) );
     SelectXovers = QQ(  randi(size(QQ,1),1) , :) 
     XoverList =  XoverList2( SelectXovers ,: ) ;
%      XoverList( randperm(size(XoverList,1) ) ,:) ;  
%      toc
   

% return
     
    fprintf('check applying %i Xovers to split the scaf routing and match the range ( %i percent) .....\n', size(XoverList,1),Acc )
    fprintf('NBaseOri~= %i  Mean: %i  Lower: %i Higher: %i .\n' ,NBaseOri ,round((NBaseOri/Npiece)) ,round((NBaseOri/Npiece)*(1-Acc/100)) ,round((NBaseOri/Npiece)*(1+Acc/100) )) ;
    AllRouting =NOriGPairBCB ;
    allGood= zeros( Npiece-1 ,1) ;
    
    
    
    
    for Nc = 1 : Npiece-1
        for k = 1:1 :  size(XoverList,1)
            OneXoverFromList= [XoverList(k,1:6) ; XoverList(k,7:12)] ;
            AllRoutingNew=GetHyperB.AddXover2OneCycle(AllRouting ,OneXoverFromList);
            
            if length(AllRoutingNew) == Nc+1  
                CheckAll =zeros( length(AllRoutingNew) ,1 ) ;
                    for cc= 1 : length(AllRoutingNew)
                    BaseRout2 = interpolateBase_ThreeCol( AllRoutingNew{cc} ); CheckAll(cc) = size(BaseRout2,1) ;
                    end                
                Good =   and( CheckAll>(NBaseOri/Npiece)*(1-Acc/100) , CheckAll<(NBaseOri/Npiece)*(1+Acc/100)) ;
                
                if sum(Good)>=Nc
                    fprintf('Cut %i : %s \n' ,Nc ,num2str(CheckAll'));              
                    allGood(Nc) = 1;
                    AllRouting= AllRoutingNew;
                    k;
                    break ;
                end
            end
        end
    end
    
        if sum(allGood)== length(allGood)  % update if every cutting is good
            GetHyperB.ScafRouting=AllRouting ;
            ConvertScafG(GetHyperB);
            ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
            ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC
            fprintf('Found the crossover satisfying the ranges, update data .....\n') ;
        else
%             figure ; ax=gca; 
%              GetHyperB.drawBCBdata2(AllRouting,157)
%              drawBCBdata2(obj,cycleList,num)
             fprintf('Cutting not good \n')
        end
    
else
    fprintf('The number of loops is not one. Nothing happens. \n')
end


 fprintf('Finish spliting Scaffold ...... \n')
end

