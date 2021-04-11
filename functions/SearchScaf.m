function SearchScaf(src,evn)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
clc;
% profile on  % last optimized, 08/26/2019 


ax= findobj(gcf,'Tag','MechScaffold3D');
axes(ax); cltab ; axes(ax);
ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;

GetHyperB.getForcedConnecteList;
GetHyperB=GetHyperB.findRT;

fH=gcf;
fbar = waitbar(0,'Please wait...'); pause(.5) ;
figure(fH);

%---------update ScafOption
nBundle=length(GetHyperB.containBundle) ;
ss_Scaffold= findobj(gcf,'Tag','ss_Scaffold');
ScafOption= ss_Scaffold.UserData.ScafOption ;
if ScafOption.choice.prevStack.scafloop.Value==1  
    opts.prevStack='scafloop' ;
else  %  PolyBundle  become meaningless
    opts.prevStack='polyT' ;
    ScafOption.choice.PolyBundle.edit.String=[] ;
end
 %------------
if ScafOption.choice.scafXover.Min.Value==1   % maximize/minimize number of scaffold Xover
    opts.scafXover='Min' ;
else  % Maximize number of scaf Xover
    opts.scafXover='Max' ;
    ScafOption.choice.NoScafXover.edit.String=[] ;     
    ScafOption.choice.minDist_btwXover.popup.Value=1 ;    
end


if isempty(ScafOption.choice.BundleNoScafLoop.edit.String);    opts.BundleNoScafLoop=[];
else
    numArr1=str2num(ScafOption.choice.BundleNoScafLoop.edit.String) ;numArr1(numArr1>nBundle)=[];
    opts.BundleNoScafLoop=numArr1 ;
    ScafOption.choice.BundleNoScafLoop.edit.String= num2str(numArr1) ;    
end
if isempty(ScafOption.choice.PolyBundle.edit.String); opts.PolyBundle=[];
else
    numArr2=str2num(ScafOption.choice.PolyBundle.edit.String) ; numArr2(numArr2>nBundle)=[];
    opts.PolyBundle=numArr2 ;
    ScafOption.choice.PolyBundle.edit.String= num2str(numArr2) ;
end
if isempty(ScafOption.choice.NoScafXover.edit.String); opts.NoScafXover=[];
else
    numArr3=str2num(ScafOption.choice.NoScafXover.edit.String) ; numArr3(numArr3>nBundle)=[];
    opts.NoScafXover=numArr3 ;
    ScafOption.choice.NoScafXover.edit.String= num2str(numArr3) ; 
end
popXX=ScafOption.choice.minDist_btwXover.popup ;
% opts.minDist_btwXover= str2num(popXX.String{popXX.Value} ) -6 ;
opts.minDist_btwXover= str2num(popXX.String{popXX.Value} ) -2 ;

GetHyperB.ScafOption= opts ;  %update to object
drawnow;
GetHyperB.CustomSkipAndInsertion = false; % whenever using routing algorithm, reset to Default
% return
%-------



%             GetHyperB.ssOption=GUIAsignSSInfo3( GetHyperB,findobj(gcf,'Tag','ss_Assembly') );   % call ssDNA asking Info GUI
fprintf(' start looking for scaffold routing \n');tic;
waitbar(.2,fbar,'start looking for scaffold routing');


scafVariety=1 ;

tol=GetHyperB.ScafOption.minDist_btwXover;
% profile on ;
switch opts.scafXover
    case 'Max'
        GetHyperB.getScafXoverAsStap  ; % get property: AllScafXover
%           GetHyperB.getScafXoverAsStap_every18  ; % get property: AllScafXover
      


        GetHyperB.findCycleList_maxScafXover([],GetHyperB.SavePremPair,2);
%         GetHyperB.findCycleList_maxScafXover([],GetHyperB.SavePremPair,2);
        
    case 'Min'
        nwhile=1;    
        while nwhile< 3  % default :5
            waitbar(0.1+nwhile*0.15 ,fbar,strcat('searching scaf routing...',num2str(nwhile) ));
%             pause(0.5)
            [~,Goodcc]=GetHyperB.findCycleList([],GetHyperB.SavePremPair,2);   %get property ScafRouting
            [ Good,whichCyl ] = CheckScafRXover( GetHyperB.ScafRouting, tol,GetHyperB.ScafOption.BundleNoScafLoop ) ;
            
            %             BundlessSpecial=[8,9,10,11];
            %             %                Bunlde2CornerZ=obj.ScafRouting((obj.ScafRouting(:,1)==2),3);Bunlde2CornerZ=unique(Bunlde2CornerZ);
            %             BunSpecialCornerZ=GetHyperB.ScafRouting(ismember(obj.ScafRouting(:,1) ,BundlessSpecial),3); BunSpecialCornerZ=unique(BunSpecialCornerZ);
            %             ExcludeRange=[70,155];   % no cross over of bundle2 locates this range,  hard coding
            %             S1= BunSpecialCornerZ<ExcludeRange(2) ;
            %             S2= BunSpecialCornerZ > ExcludeRange(1);
            %             SS= and(S1,S2);            
            nwhile=nwhile+1;
            fprintf('while loop iteration for scafold is  %d, %i \n',nwhile,size(whichCyl,1))            
            if Good==1  || size(whichCyl,1)==2  || Goodcc==1 % && sum(SS)==0
                break;
            end
        end
        Good
        whichCyl
end
% profile viewer ;

GetHyperB.Scaf_fromCadDOM= GetHyperB.ScafRouting ;
toc;

ConvertScafG(GetHyperB);
ConvertScafSQ(GetHyperB);      % Get properties:ScafdigitSQ
ConvertScafHC(GetHyperB);      % Get properties:ScafdigitHC

fprintf(' found scaffold routing Again\n')
waitbar(.9,fbar,'found scaffold routing Again');

axes(ax);
hold on; axis  equal;
% return
% hSurf= GetHyperB.plotScafR_cylindermodel ;

 hSurf= GetHyperB.plotScafR_cylindermodelMulti ;

%             ForLegend=surface(nan, nan,'Tag','test','EdgeColor','interp','FaceColor','none','Tag','DefaultScafH');
            % ForLegend=surface(nan, nan, 'Linestyle', 'none', 'Marker', 'none' );
%                         ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');

            hLg= legend(hSurf{1},'Click me for instructions','Location','northwest' ) ; hLg.String={' \rmScaffold (index \bf\color[rgb]{0.24,0.15,0.66}Start 5'' \color{black}to \color[rgb]{0.78,0.78,0.18}End 3''\rm\color{black})'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
%             ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_Scaf( src,evn,ax );
            hLg.Title.String='Click me for instructions' ;
            hLg.Units='normalized'; %hLg.AutoUpdate ='off';
            hLg.Position=[0.0063 0.9528 0.1569 0.0387];

close(fbar) ;
% profile viewer

% \color[rgb]{0,0.5,0}

end

