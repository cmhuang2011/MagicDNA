


if 0==1
ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;  % select the program with the design as routing template
GetHyperB1= ss_Assembly.UserData.HyperBundle ;

ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;  % select the program with the design want to substitute
GetHyperB2= ss_Assembly.UserData.HyperBundle ;
end

%---------------scaffold 

GetHyperB2.ScafRouting{1} =GetHyperB1.ScafRouting{1}  ;




%---------------staple 

[BundlesIndex_M, TF_CrossBundles_M]=GetHyperB1.ShowStapleInBundles;
[a_M,b_M ]=ismember(BundlesIndex_M, [1:3] ) ;   % enter bundle indexes 
CopyStap = GetHyperB1.StapList3(a_M) ;
 
[BundlesIndex, TF_CrossBundles]=GetHyperB2.ShowStapleInBundles;
[a,b ]=ismember(BundlesIndex, [1:3] ) ;
OriStap3_S=GetHyperB2.StapList3 ;
 
OriStap3_S=OriStap3_S(~a) ;
New =[OriStap3_S ;CopyStap] ;
GetHyperB2.StapList3=New ;


%------------
ConvertScafG(GetHyperB2);
ConvertScafSQ(GetHyperB2);      % Get properties:ScafdigitSQ
ConvertScafHC(GetHyperB2);      % Get properties:ScafdigitHC


GetHyperB2.ConvertStap('Square');    % Get properties: DigitStapSQ,   HeadOfStep
GetHyperB2.ConvertStap('Honeycomb');    % Get properties: DigitStapHC


