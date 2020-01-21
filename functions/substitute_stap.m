function substitute_stap( SourceHB, TargetHB, SourceBundle, TargetBundle )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

S_RT=SourceHB.RelateTable ;
T_RT=TargetHB.RelateTable ;

S_Cylinders = S_RT(ismember(S_RT(:,1),SourceBundle) ,:) ;

T_Cylinders = T_RT(ismember(T_RT(:,1),TargetBundle) ,:) ;
CombineTable = [S_Cylinders , T_Cylinders] ;   % if cylinder number not match, expect to have error

C5Mapping = CombineTable(:,[1 5  8 12]) ;

[BundlesIndex_T, ~]=TargetHB.ShowStapleInBundles ;
TargetStapList3 = TargetHB.StapList3 ;
TargetStap_toRemove = ismember( BundlesIndex_T, C5Mapping(:,3))   ;
TargetStap_toRemain = ~TargetStap_toRemove   ;


SourceStapList3 = SourceHB.StapList3 ;
[BundlesIndex_S, ~]=SourceHB.ShowStapleInBundles ; 
SourceStap_toAdd_OriCyl =  ismember( BundlesIndex_S, C5Mapping(:,1))   ;   % ind


%-+-------mapping cylinder index 
 OldStap3_toAdd = SourceStapList3(SourceStap_toAdd_OriCyl) ;
for k =1:length(OldStap3_toAdd)
 Old = OldStap3_toAdd{k} ;
 
 [~,bb] =ismember(Old(:,1), C5Mapping(:,2) ) ;
 NewInds =  C5Mapping(bb,4) ;
 NewStap = [NewInds, Old(:,2)] ;
 
 OldStap3_toAdd{k}=NewStap ;
end

TargetHB.StapList3 =    [TargetStapList3(TargetStap_toRemain);OldStap3_toAdd ] ;


% sdfsf=2

% SourceStapList3{SourceStap_toAdd}


TargetHB.ConvertStap('Square');    % Get properties: DigitStapSQ,   HeadOfStep
TargetHB.ConvertStap('Honeycomb');    % Get properties: DigitStapHC




% [BundlesIndex_M, TF_CrossBundles_M]=GetHyperB_S.ShowStapleInBundles ;


% sdfs=3

end

