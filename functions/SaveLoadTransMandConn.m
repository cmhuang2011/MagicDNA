function Save=SaveLoadTransMandConn(IsSave  )
%%% CM use, recover transformation matrix and connectivity

%   Detailed explanation goes here

getHBundleHandle ;   %  get handle GetHyperB;
Save=[];
if IsSave==1
    Save.SaveTransM = cell(size(GetHyperB.containBundle)) ;
    for k=1: length(GetHyperB.containBundle )       
        Save.SaveTransM{k} = GetHyperB.containBundle{k}.TransformMatrix2 ;
    end
    Save.SaveAdjM =GetHyperB.BundleAdjM ;
    save('Temp.mat','Save') ;    
    disp('The adjacency matrix between bundles :') ;
    disp(GetHyperB.BundleAdjM) ; 
    disp('have been saved in Temp.mat.') ;
    
    G=graph(GetHyperB.BundleAdjM) ; 
    disp('Adjacency list: ') ; 
    disp(G.Edges )
%     A = adjacency(G)
    disp('The transformation matrixs are also saved for revising geometries in the table in STEP.') ;
    
    
    
else  % load     
    load('Temp.mat','Save')  ;    
end
end

