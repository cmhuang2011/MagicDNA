function  DelAddBundle( src,evn,OldHyperB,fH,patchH,popupH )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rotate3d off ;    % important, prevent users add/delete bundle with rotate3D on


% NewFig=figure;

NewFig = dialog('units','pixels',...
    'position',[300 300 400 500],...
    'menubar','none',...
    'name','Add or Remove Bundles',...
    'numbertitle','off',...
    'resize','off');

movegui(NewFig,'center')
N=length(OldHyperB.containBundle);


Mat=[ (1:N)'  ,zeros(N,1)]  ;
Datat=mat2cell(Mat, ones(1,N ) , [1,1]  );
for rev=1:size(Datat,1)
    Datat{rev,1}= num2str(Datat{rev,1});
    Datat{rev,2}=true;
end
RedSelectBun = popupH.Value; 
if length(RedSelectBun)>1
    for rev=1:length(RedSelectBun)
%      Datat{RedSelectBun(rev),2}=false;
    end
%     Datat{rev,2}=true; 
end


Datat{end+1,1}='New Bundle/Mech.' ;
Datat{end,2}=' ' ;


t = uitable(NewFig,'Units','normalized','Position',[0.05 0.2 0.9 0.75],....
    'ColumnWidth','auto','ColumnFormat',({[] [] }),...
    'ColumnEditable', true,'Data',Datat,....
    'ColumnName',{'Old Bundle'; 'Keep/Delete' });

btn_AddBundle = uicontrol('Style', 'pushbutton', 'String', '	Add New Bundle/Mech. ','Unit','normalized', 'Position', [0.05 0.05 0.4 0.1] );
btn_ExportNewHB = uicontrol('Style', 'pushbutton', 'String', '	OK ','Unit','normalized', 'Position', [0.5 0.05 0.4 0.1] );


NewFig.UserData.Extra=[];
btn_AddBundle.Callback=@(src,evn)AddBun(src,evn,t) ;

btn_ExportNewHB.Callback=@(src,evn)ExportHB(src,evn,t,OldHyperB,fH) ;


% foo = uitable('Data', {false;true;true;false;false}, 'ColumnEdit', true);


% EraseB=3;
% OldBunds=1:length(OldHyperB.containBundle) ;
% RemainBundel=setdiff(OldBunds,EraseB) ;
%
% NewAdjM= fH.UserData.BundleAdj(RemainBundel,RemainBundel);
%
%
% InportHyperBundle=[];
% InportHyperBundle=hyperbundle(1,OldHyperB.containBundle{RemainBundel(1) } );
%
% for k=2:length(RemainBundel)
%     AddBundle=  OldHyperB.containBundle{RemainBundel(k) } ;
%     InportHyperBundle=AddBundles(InportHyperBundle,AddBundle);
% end
%
%
% prefercase = AskPrePairPreference();
%
% trial=0;
% while trial<=1000
% [Doable,CellPairList,CellUnpair]=checkpairableH(InportHyperBundle,prefercase.pair);
%  if  (Doable==1   ||  trial>=300)
%        break;
%  end
%  trial=trial+1;
% end
%
% PremPair.Doable=Doable;
% PremPair.CellPairList=CellPairList;
% PremPair.CellUnpair=CellUnpair;
%
%
% InportHyperBundle=InportHyperBundle.CreateTFWindow(PremPair,NewAdjM) ;
%
end

function AddBun(src,evn,t)
fH2=gcf;
PartFolder=[ pwd '\parts\'];
%     filename = uigetfile_with_preview( '*.mat')
[FileName,PathName] = uigetfile('*.mat','DialogTitle',PartFolder);
SSS=open(strcat(PathName,FileName));
%           uiopen('load')
if  isfield(SSS,'Psave')   % user select a  part
    TempPart=SSS.Psave.GUISavePart ;
    fH2.UserData.Extra{end+1}=TempPart;
    fH2.UserData.Extra   ;
    t.Data{end+1,1}=FileName(1:end-4) ;
    t.Data{end,2}=true;
    
elseif isfield(SSS,'S') % user select a  mechanism(one of multiple parts)
    %      sdffs=3
    IndeMap = [-1,-1];
    for k=1: length(SSS.S.containBundle)
        TempPart=  SSS.S.containBundle{k} ;
        fH2.UserData.Extra{end+1}=TempPart;
        fH2.UserData.Extra   ;
        IndeMap=union(IndeMap,[size(t.Data,1) ,k] ,'rows') ;
        t.Data{end+1,1}=strcat(FileName(1:end-4) ,'-' ,num2str(k) ) ;
        t.Data{end,2}=true;
    end
    IndeMap=setdiff(IndeMap, [-1,-1],'rows') ;
    [u,v] =find( SSS.S.BundleAdjM) ;  % change to list
    Q= u<v;
    u=u(Q) ; v=v(Q) ; W= zeros(size(v)) ;
    for k2= 1:length(W)
        W(k2) =  SSS.S.BundleAdjM(u(k2), v(k2)) ;
    end
    [a1,b1]=ismember(u,IndeMap(:,2)) ;
    [a2,b2]=ismember(v,IndeMap(:,2)) ;
    List= [ IndeMap(b1,1), IndeMap(b2,1),W] ;  % correct index
    if isfield(fH2.UserData,'AddAdjL' )
        fH2.UserData.AddAdjL = [ fH2.UserData.AddAdjL ;List];
    else
        fH2.UserData.AddAdjL =List    ;
    end
    
            dsdf=3
end


end

function ExportHB(src,evn,t,OldHyperB,fH)
fH2=gcf ;
% tic
Nold=   length(OldHyperB.containBundle);
OldBunds=1:length(OldHyperB.containBundle) ;

EraseB=[];
for k=1:length(OldBunds)
    if t.Data{OldBunds(k),2}==0
        EraseB=union(EraseB,k) ;
    end
end

RemainBundel=setdiff(OldBunds,EraseB) ;
%
NewAdjM= fH.UserData.HyperBundle.BundleAdjM(RemainBundel,RemainBundel);
%
%
InportHyperBundle=[];


Z1Arr=OldHyperB.containBundle{RemainBundel(1) }.Zbase1 ;
Z2Arr=OldHyperB.containBundle{RemainBundel(1) }.Zbase2 ;
InplaneXY = OldHyperB.containBundle{RemainBundel(1) }.CylInplanePosition ;
    if strcmp(OldHyperB.containBundle{RemainBundel(1) }.type,'SQ')
        AddBundle=  BundleCylinderSQ(1,[],Z1Arr,Z2Arr,InplaneXY) ;
    else
        AddBundle=  BundleCylinderHC(1,[],Z1Arr,Z2Arr,InplaneXY) ;
    end
AddBundle.TransformMatrix2=OldHyperB.containBundle{RemainBundel(1) }.TransformMatrix2;

InportHyperBundle=hyperbundle(1,AddBundle );
% InportHyperBundle=hyperbundle(1,OldHyperB.containBundle{RemainBundel(1) } );

for k=2:length(RemainBundel)
%     AddBundle=  OldHyperB.containBundle{RemainBundel(k) } ;
    
    Z1Arr=OldHyperB.containBundle{RemainBundel(k) }.Zbase1 ;
    Z2Arr=OldHyperB.containBundle{RemainBundel(k) }.Zbase2 ;
    InplaneXY = OldHyperB.containBundle{RemainBundel(k) }.CylInplanePosition ;
    if strcmp(OldHyperB.containBundle{RemainBundel(k) }.type,'SQ')
        AddBundle=  BundleCylinderSQ(1,[],Z1Arr,Z2Arr,InplaneXY) ;
    else
        AddBundle=  BundleCylinderHC(1,[],Z1Arr,Z2Arr,InplaneXY) ;
    end
    AddBundle.TransformMatrix2=OldHyperB.containBundle{RemainBundel(k) }.TransformMatrix2;
    
    
    InportHyperBundle=AddBundles(InportHyperBundle,AddBundle);
end

for kk=Nold+2:size( t.Data,1)
    if t.Data{kk,2}==1
        AddBundle=  fH2.UserData.Extra{kk -Nold-1}       ;
        InportHyperBundle=AddBundles(InportHyperBundle,AddBundle);   NH= length(InportHyperBundle.containBundle) ;
        NewAdjM(end+1,:)=0;
        NewAdjM(:,end+1)=0;
        %     InportHyperBundle.containBundle{NH}.
    end
end
% gcf
EraseB2=[];    % if users add mechanism, also import
for k=max(OldBunds)+2:size(t.Data,1)
    if t.Data{k,2}==0
        EraseB2=union(EraseB2,k-1) ;
    end
end
RemainBundel2=setdiff(1:size(t.Data,1)-1,union(EraseB,EraseB2) ) ;
if isfield(fH2.UserData,'AddAdjL')
    for k=1:size(fH2.UserData.AddAdjL ,1)
        if ismember(fH2.UserData.AddAdjL(k,1) , RemainBundel2) && ismember(fH2.UserData.AddAdjL(k,2) , RemainBundel2)
            NewInd1 = fH2.UserData.AddAdjL(k,1)==RemainBundel2  ;
            NewInd2 = fH2.UserData.AddAdjL(k,2)==RemainBundel2  ;
            NewAdjM(  NewInd1,NewInd2)=fH2.UserData.AddAdjL(k,3) ;
            NewAdjM(  NewInd2,NewInd1)=fH2.UserData.AddAdjL(k,3) ;   %-----
            
        end
    end
end
close(fH2); % Closing the table here will erase the sound. 0719

% NewAdjM=zeros(size(NewAdjM)) ; % hard
%
prefercase = AskPrePairPreference ;
trial=0;
while trial<=1000
    trial=trial+1;
    [Doable,CellPairList,CellUnpair]=checkpairableH(InportHyperBundle,prefercase.pair);
    if  (Doable==1   ||  trial>=500)
        break;
    end
end
% [Doable ,trial]
if ~isempty(OldHyperB.SavePremPairManual) && isempty(setdiff(RemainBundel,RemainBundel2))  && isempty(setdiff(RemainBundel2,RemainBundel))  % no adding or remove bundles
  for k= 1 : length(RemainBundel2)    
      CellPairList{k} = OldHyperB.SavePremPairManual.CellPairList{RemainBundel2(k)} ;
      CellUnpair{k} = [];
  end
%     sdfsf=3
end


CellUnpair;
BundleGroupWithSameCylPosition = InportHyperBundle.InplanePairingMatch  
for i= 1 :max(BundleGroupWithSameCylPosition)
    IncludeBundles= find(BundleGroupWithSameCylPosition == i) ;
    for k=2: length(IncludeBundles)
        CellPairList{IncludeBundles(k)} =CellPairList{IncludeBundles(1)} ;
    end
end
fprintf('Force bundles with same cylinderPosition using same pairing. \n')
% sdf=3
%-----------------
% fw=msgbox('Using hard code for pairing !! ');
% waitfor(fw);
% fprintf('Using hard code for pairing !! \n')
% % CellPairList{1}=[1 7;2 8; 3 4; 5 9; 6 10; 11 13 ; 15 17 ;12 14;16 18; 19 23; 20 24; 25 26; 21 27 ; 22 28]  ;
% % CellPairList{7}=[1 2 ; 5 9;6 10;13 14 ;4  8 ;11 12; 3 7]  ;
% % GroupPair=setdiff(2:12,14);% hard code
% GroupPair=2:12;% hard code
% for k=1:length(GroupPair)    % hard code
%    CellPairList{GroupPair(k)}=CellPairList{1} ;
% end
%-------------------
% %
% fprintf('Using hard code for pairing !! \n')
% GroupPair2=4:6 ;  % hard code
% for k=1:length(GroupPair2)
%  CellPairList{GroupPair2(k)}=[2,1 ;4,6;5,3] ;
% end
%---------------

fprintf(' Doable =%i , trial= %i \n',Doable,trial) ;
PremPair.Doable=Doable;
PremPair.CellPairList=CellPairList;
PremPair.CellUnpair=CellUnpair;

axes(findobj(fH,'Tag','AssemblyMain')) ;
cltab;
% toc
InportHyperBundle=InportHyperBundle.CreateTFWindow(PremPair,NewAdjM,[],fH) ;
% toc
% close(fH2);  % Closing the table here will make some sound. 0719


end

% end

