GetHyperB2.ScafRouting{1} =GetHyperB1.ScafRouting{1}  ;


[BundlesIndex_M, TF_CrossBundles_M]=GetHyperB1.ShowStapleInBundles;
[a_M,b_M ]=ismember(BundlesIndex_M, [1:3] ) ;
CopyStap = GetHyperB1.StapList3(a_M) ;

[BundlesIndex, TF_CrossBundles]=GetHyperB2.ShowStapleInBundles;
[a,b ]=ismember(BundlesIndex, [1:3] ) ;
GetHyperB2.StapList3(~a)

OriStap3_S=OriStap3_S(~a) ;
New =[OriStap3_S ;CopyStap] ;
GetHyperB2.StapList3=New ;


