function ReOrderBundles( OldHyperB,ftab )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

rotate3d off ;    % important, prevent users add/delete bundle with rotate3D on


% NewFig=figure;

NewFig = dialog('units','pixels',...
    'position',[300 300 400 500],...
    'menubar','none',...
    'name','Export a new assemblt with reordered bundles',...
    'numbertitle','off',...
    'resize','off');


movegui(NewFig,'center')
N=length(OldHyperB.containBundle) ;

Mat=[ (1:N)'  ,zeros(N,1)]  ;
Datat=mat2cell(Mat, ones(1,N ) , [1,1]  );
for rev=1:size(Datat,1)
    Datat{rev,1}= num2str(Datat{rev,1});
    Datat{rev,2}=rev;
end


% 
% Datat{end+1,1}='New Bundle/Mech.' ;
% Datat{end,2}=' ' ;
% 
% Datat{end+1,1}='New Bundle/Mech.' ;
% Datat{end,2}=' ' ;

str = cellfun(@num2str,num2cell(1:N),'uniformoutput',0);

% ,'ColumnFormat',({[] [] strHC str})
t = uitable(NewFig,'Units','normalized','Position',[0.05 0.2 0.9 0.75],....
    'ColumnWidth','auto','ColumnFormat',({[] str }),...
    'ColumnEditable', true,'Data',Datat,....
    'ColumnName',{'Old Bundle index'; 'New Bundle index' });

Export_btn = uicontrol('Style', 'pushbutton', 'String', '	Export  ','Unit','normalized', 'Position', [0.05 0.05 0.4 0.1] );

Export_btn.Callback=@(src,evn)ExportNewAsOrder(src,evn,t,OldHyperB,ftab,NewFig) ;
% sdsf=3

end

function ExportNewAsOrder(src,evn,t,OldHyperB,ftab,NewFig) 

OldOrders=  str2num(cell2mat( t.Data(:,1) )) ;
NewOrders=  cell2mat( t.Data(:,2) ) ;

if length(intersect(NewOrders,OldOrders))~= length(OldOrders)

%     union(NewOrders,OldOrders)
    f = msgbox('Invalid input. Please assign an array with all numbers ');
%     NotGood =1
  
else
    %      sdf=3
    close(NewFig); % Closing the table here will erase the sound. 0719

    OldAdjM=  OldHyperB.BundleAdjM ;
    [a,b]= ismember(OldOrders,NewOrders) ;
    NewAdjM=OldAdjM(b , b ) ;
    %      find(NewOrders )
    InportHyperBundle=[];
    InportHyperBundle=hyperbundle(1,OldHyperB.containBundle{find(NewOrders==1) } );
    
    for k=2: length(NewOrders)
        AddBundle=  OldHyperB.containBundle{find(NewOrders==k) } ;
        InportHyperBundle=AddBundles(InportHyperBundle,AddBundle);   
              
    end
    
    prefercase = AskPrePairPreference ;
    trial=0;
    while trial<=1000
        trial=trial+1;
        [Doable,CellPairList,CellUnpair]=checkpairableH(InportHyperBundle,prefercase.pair);
        if  (Doable==1   ||  trial>=500)
            break;
        end
    end
    
    BundleGroupWithSameCylPosition = InportHyperBundle.InplanePairingMatch
    for i= 1 :max(BundleGroupWithSameCylPosition)
        IncludeBundles= find(BundleGroupWithSameCylPosition == i) ;
        for k=2: length(IncludeBundles)
            CellPairList{IncludeBundles(k)} =CellPairList{IncludeBundles(1)} ;
        end
    end
    fprintf('Force bundles with same cylinderPosition using same pairing. \n')
    
    fprintf(' Doable =%i , trial= %i \n',Doable,trial) ;
    PremPair.Doable=Doable;
    PremPair.CellPairList=CellPairList;
    PremPair.CellUnpair=CellUnpair;
    
    axes(findobj(ftab,'Tag','AssemblyMain')) ;
    cltab;
    % toc
    InportHyperBundle=InportHyperBundle.CreateTFWindow(PremPair,NewAdjM,[],ftab) ;
    
    %
    % for k=2:length(RemainBundel)
%     AddBundle=  OldHyperB.containBundle{RemainBundel(k) } ;
%     InportHyperBundle=AddBundles(InportHyperBundle,AddBundle);
% end 
   
    
end

end


