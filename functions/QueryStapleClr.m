function  Result=QueryStapleClr( GetHyperB, MagicDNAfH)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Result =[];
% ss_Assembly= findobj(0,'Tag','ss_Assembly') ;
% GetHyperB= ss_Assembly.UserData.HyperBundle ;

rotate3d off ;    % important, prevent users add/delete bundle with rotate3D on
DefaultForAll = 6 ;

global minDist_StapleXoverFromTwoSide

if ~isempty(minDist_StapleXoverFromTwoSide)
    if size(minDist_StapleXoverFromTwoSide,1) == length(GetHyperB.containBundle)
    DefaultArray = minDist_StapleXoverFromTwoSide ;
    else
    DefaultArray =  DefaultForAll*ones(length(GetHyperB.containBundle),2) ;   
    end
else
    DefaultArray =  DefaultForAll*ones(length(GetHyperB.containBundle),2) ;   
end

% return
% NewFig=figure;

NewFig = dialog('units','pixels',...
    'position',[300 300 400 500],...
    'menubar','none',...
    'name','Specify individual min bps to exclude staple Xovers from bundle ends',...
    'numbertitle','off',...
    'resize','off');

movegui(NewFig,'center')
N=length(GetHyperB.containBundle);
% N = 12 ;

Mat=[ (1:N)'  ,zeros(N,1)]  ;
Datat=mat2cell(Mat, ones(1,N ) , [1,1]  );
for rev=1:size(Datat,1)
    Datat{rev,1}= num2str(DefaultArray(rev,1));
    Datat{rev,2}= num2str(DefaultArray(rev,2));
%     Datat{rev,1}= num2str(DefaultForAll);
%     Datat{rev,2}= num2str(DefaultForAll);
    
end


Datat{end+1,1}=num2str(DefaultForAll);
Datat{end,2}=num2str(DefaultForAll);
str = cellfun(@num2str,num2cell([2:10 , 15]),'uniformoutput',0);

t = uitable(NewFig,'Units','normalized','Position',[0.05 0.2 0.9 0.75],....
    'ColumnWidth','auto','ColumnFormat',({str str }),...
    'ColumnEditable', true,'Data',Datat,....
    'ColumnName',{'Z1 site'; 'Z2 site' });

BaColor =  repmat(t.BackgroundColor, 2+floor(N/2) ,1) ;
BaColor=BaColor(1:N+1,:) ;
BaColor(end,:) = [0.3010, 0.7450, 0.9330] ;

t.BackgroundColor=BaColor ;


btn_ShowInd = uicontrol(NewFig,'Style', 'pushbutton', 'String', 'Show bundle indexes','Unit','normalized', 'Position', [0.05 0.05 0.4 0.1] );
btn_OK = uicontrol(NewFig,'Style', 'pushbutton', 'String', '	OK ','Unit','normalized', 'Position', [0.5 0.05 0.4 0.1] );


t.CellEditCallback=  @(src,evn)TableECK(src,evn);



% NewFig.UserData.Extra=[];
btn_ShowInd.Callback=@(src,evn)ShowBundleInd(src,evn ,MagicDNAfH) ;
btn_OK.Callback=@(src,evn)Export(src,evn,t,NewFig) ;

uiwait(NewFig)
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

function Export(src,evn,t,NewFig)
% Result= t.Data; 

QQ = zeros(size(t.Data,1)-1 , 2) ;
for k =1:size(QQ,1)
    QQ(k,1) = str2double(t.Data{k,1}) ;
    QQ(k,2) = str2double(t.Data{k,2}) ;   
end
Result=QQ ;
close(NewFig); % Closing the table here will erase the sound. 0719

end

end

function TableECK(src,evn)
if strcmp(evn.EventName,'CellEdit')
    if evn.Indices(1)==size(src.Data,1) %  change last row
        switch  evn.Indices(2)
            case 1
                for k=1:size(src.Data,1)-1
                    src.Data{k,1}= src.Data{end,1} ;
                end
            case 2
                for k=1:size(src.Data,1)-1
                    src.Data{k,2}= src.Data{end,2} ;
                end
        end
    end
end
end



function ShowBundleInd(src,evn,MagicDNAfH)
ax=findobj(MagicDNAfH,'Tag','AssemblyMain');
% axes(ax) ;

% SaveAxes_two ;

ax_old = ax;

f_new = figure; 




ax_new = copyobj(ax,f_new)  ;
set(ax_new,'Position','default');
ax_new.Tag='';
set(gca,'FontSize',20)

% fprintf('Use ex: print(''-r100'',f_new,''Image3_r100'',''-djpeg'') \n' ) ;
% fprintf('Use ex: set(gcf,''KeyPressFcn'',@(src,evn)ChangeAxesLimit(src,evn) ) \n' ) ;

if isfield(ax_old.UserData ,'hlink')
ax_old_hidden = ax_old.UserData.hlink.Targets(2) ;
ax_new2 =  axes('Position',ax_new.Position,'Visible','off','hittest','off'); % Invisible axes
% ax_new2 =  axes('Parent', f_new,'Position',ax_new.Position,'Visible','off','hittest','off'); % Invisible axes
% ax_new2 =axes('Position',ax_new.Position,'Visible','off','hittest','off') ;
h_text = findobj(ax_old_hidden,'Type','Text')  ;
h_text_new = copyobj(h_text,ax_new2)  ;

hlink =linkprop([ax_new ax_new2],{ 'XLim' 'YLim' 'ZLim' ,'PlotBoxAspectRatio','View' }); % The axes should stay aligned
ax_new.UserData.hlink= hlink ;
% axes(ax_new) ;
end



end


