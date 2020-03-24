%% This is the script for designing the 3 scaffold plates with overhangs.
% The goal is to use overhang locations to show specific text. In this case
% it is 'O-S-U' three characters.
% Three plates are connected by staple overhangs. and we want to have
% common staples except staples with overhans.

% asseumption: three identical cylinder model(double-stranded area)


switch stepAssign
    %%  Step 1
    % make three scaffold share the same routing.
    case 1
        % run this section by copy-paste the following to command line 
%                    stepAssign=1 ;
%                    Ex_DecoratePlatesWithOverhangs ;

        
        ss_Assembly= findobj(0,'Tag','ss_Assembly') ;
        GetHyperB= ss_Assembly.UserData.HyperBundle ;   % get hyperbundle handle
        
        
        %check original random routing :
        GetHyperB.ScafRouting
        
        %BCB representation: simply change bundle indexes in the other two ;
        GetHyperB.ScafRouting{2}=GetHyperB.ScafRouting{1} ; GetHyperB.ScafRouting{2}(:,1)=2;
        GetHyperB.ScafRouting{3}=GetHyperB.ScafRouting{1} ;  GetHyperB.ScafRouting{3}(:,1)=3;
        
        
        %----convert to digital format for cadnano
        ConvertScafG(GetHyperB);
        ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
        ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC
        
        %---------
        % Use UI to visualize the three scaffold routing
        

        
    case 2
        % Step 2
        % The staple overhang algorithm basically has two parts: 1. leave nicks on
        % the positions in  GetHyperB.FindStap. 2. extand overhang in
        % GetHyperB.extendOverhang
        
        % ---------
        % In this case, we want leave nicks on three top surfaces and side
        % connections in the first parts.
        % in GetHyperB.FindStap
        % 
        InitialStappEnds   ; % C3 notation   % input
        
        BCB_stapEnds = zeros(size(InitialStappEnds,1) ,3) ;
        BCB_stapEnds(:,3) = InitialStappEnds(:,2) ;
        BCB_stapEnds(:,1:2) =  obj.RelateTable( InitialStappEnds(:,1),1:2)  ;
        
        sum(BCB_stapEnds(:,1)==1) ;
        sum(BCB_stapEnds(:,1)==2) ;
        sum(BCB_stapEnds(:,1)==3) ;
        
        Common_CB = unique(BCB_stapEnds(:,2:3),'rows') ;
        BCB1= [ones(size(Common_CB,1) ,1) , Common_CB] ;
        BCB2= [2*ones(size(Common_CB,1) ,1) , Common_CB] ;
        BCB3= [3*ones(size(Common_CB,1) ,1) , Common_CB] ;
        New_BCB_all =[BCB1;BCB2 ;BCB3] ;
        % -----convert back to C3
        [aa,bb] =ismember(New_BCB_all(:,1:2) , obj.RelateTable(:,1:2),'rows') ;
        NewC3 = [bb,New_BCB_all(:,3)] ;
        
        InitialStappEnds=  NewC3 ; % output        
    case 3
       % continue to apply crossovers and breaking staples
       % later copy bundle 2 staple routing to bundle 1 and 3, before
       % extanding overhang
       % 
       % in SearchStap
       
       Staple5pInBundles =  GetHyperB.ShowStapleInBundles     ;   % obj.StapList3 ; % C5 Rep
       OriAll_StapList3 = GetHyperB.StapList3 ; 
       Ind = Staple5pInBundles==2 ;
       StapListB2 = OriAll_StapList3(Ind) ; % C5 representation 
       StapListB1= StapListB2 ;
       for k =1:length(StapListB1)
           QQ= StapListB1{k} ;
           [aa,bb] =ismember(QQ(:,1) , GetHyperB.RelateTable(:,5),'rows') ;
           BCinBundle2 =GetHyperB.RelateTable(bb,1:2) ;
           BCinBundle1 =BCinBundle2 ; BCinBundle1(:,1) =1 ;
           
           [aa2,bb2] =ismember(BCinBundle1 , GetHyperB.RelateTable(:,1:2),'rows') ;
           StapListB1{k}=[ GetHyperB.RelateTable(bb2,5), QQ(:,2) ];
       end
       StapListB3= StapListB2 ;
       for k =1:length(StapListB3)
           QQ= StapListB3{k} ;
           [aa,bb] =ismember(QQ(:,1) , GetHyperB.RelateTable(:,5),'rows') ;
           BCinBundle2 =GetHyperB.RelateTable(bb,1:2) ;
           BCinBundle3 =BCinBundle2 ; BCinBundle3(:,1) =3 ;
           
           [aa2,bb2] =ismember(BCinBundle3 , GetHyperB.RelateTable(:,1:2),'rows') ;
           StapListB3{k}=[ GetHyperB.RelateTable(bb2,5), QQ(:,2) ];
       end       
       GetHyperB.StapList3= [StapListB1 ; StapListB2;StapListB3];
    case 4
        % let the extanding overhang continue without any hard coding since
        % we have left nicks (more than required) on those positions. The
        % table in overhang design should guide the extension on sides as
        % connections and on surfaces as 'O' 'S' 'U"
       
        % export staple sequences and check number of unique staples
        size(unique( CheckRepeatSeq.Seq))
        size( CheckRepeatSeq.Seq) 
        
    case 5
        % For oxDNA simulation, we don't want the closing strands on the
        % three surfaces but keep the closing strands for connection. 
        % Just keep first 12 pairs closing strands .
       Inds =1:12 ;
       ClosingStrandL    =ClosingStrandL(Inds) ;
       
       % hard code for  the following in oxDNAInitail.ExportOxDNA;
       % BM2=[ cell2mat(scafBundleRout); cell2mat(stapBundleRout);cell2mat(CSBundleRout(1:12) )] ;
       
end
