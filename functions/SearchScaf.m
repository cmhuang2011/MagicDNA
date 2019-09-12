function SearchScaf(src,evn)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
clc;
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
opts.minDist_btwXover= str2num(popXX.String{popXX.Value} ) -2 ;
GetHyperB.ScafOption= opts ;  %update to object
drawnow;
% return
%-------



%             GetHyperB.ssOption=GUIAsignSSInfo3( GetHyperB,findobj(gcf,'Tag','ss_Assembly') );   % call ssDNA asking Info GUI
fprintf(' start looking for scaffold routing \n');tic;
waitbar(.2,fbar,'start looking for scaffold routing');


scafVariety=1 ;

tol=GetHyperB.ScafOption.minDist_btwXover;

switch opts.scafXover
    case 'Max'
        GetHyperB.getScafXoverAsStap  ; % get property: AllScafXover
%           GetHyperB.getScafXoverAsStap_every18  ; % get property: AllScafXover
      
        %----------detect non-boundary cylinder, take out a portion of scaf Xover for inner cylinder
%         BundleCylinders= [-1, -1] ;
%         for Bi=1 :length(GetHyperB.containBundle)
%             ExternalXoversOnIt= GetHyperB.containBundle{Bi}.ExternalXoverAsFB(:,1) ;ExternalXoversOnIt(ExternalXoversOnIt<0) =[];
%             ExternalXoversOnIt=unique(ExternalXoversOnIt);
%             BundleCylinders= union(BundleCylinders, [Bi*ones(size(ExternalXoversOnIt)),ExternalXoversOnIt],'rows' ) ;
%         end
%         BundleCylinders=setdiff(BundleCylinders,[-1,-1] ,'rows') ;
%         IndArrange1 = 1:2:size(GetHyperB.AllScafXover,1) ;
%         IndArrange2 = 2:2:size(GetHyperB.AllScafXover,1) ;
%         AllScafXoverInOneRow =[  GetHyperB.AllScafXover(IndArrange1,:),   GetHyperB.AllScafXover(IndArrange2,:)] ;
%         fprintf('total number of scaf Xover before algorithm = %i \n' ,  size(AllScafXoverInOneRow,1));
%         
%         IsInternalXover1 = ~ismember(AllScafXoverInOneRow(:,1:2), BundleCylinders,'rows') ;
%         IsInternalXover2 = ~ismember(AllScafXoverInOneRow(:,4:5), BundleCylinders,'rows') ;
%         IsInternalXover3 = ~ismember(AllScafXoverInOneRow(:,7:8), BundleCylinders,'rows') ;
%         IsInternalXover4 = ~ismember(AllScafXoverInOneRow(:,10:11), BundleCylinders,'rows') ;
%         %-definition of internal : Both side are internal cylinders(and)
%         IndInternal= and(and(and(IsInternalXover1,IsInternalXover2) ,IsInternalXover3) ,IsInternalXover4) ;
%         fprintf('total number of Internal scaf Xover= %i \n' ,  sum(IndInternal));
%         %-----------
%         IntXover =AllScafXoverInOneRow(IndInternal,:)   ;
%         RatioOut =0.3 ; fprintf('Ratio of removed internal Xover %3.2f \n',RatioOut ) ;
%         while 1
%             Inds =rand( sum(IndInternal),1) <RatioOut ;
%             if sum( Inds(1:end-1)+Inds(2:end)==2)==0 && abs(sum(Inds)/length(Inds)-RatioOut)<0.05
%                 break % non continous
%             end
%         end
%         fprintf('following is removed Xover  \n' )  ;
%         IntXover(Inds ,:)        
%         IntXover=IntXover(~Inds ,:) ;  % remove random pick ones        
%         RestXover = [IntXover ;AllScafXoverInOneRow(~IndInternal,:)  ] ;
%         fprintf('total number of rest scaf Xover  = %i \n' ,  size(RestXover,1));
%         
%         RestXover2=[RestXover(:,1:6) ;RestXover(:,7:12) ];
%         Indd1= (1:1:size(RestXover,1));Indd1=[Indd1;Indd1+max(Indd1)];
%         Indd1=reshape(Indd1,[ 2*size(Indd1,2),1]) ;
% %         GetHyperB.AllScafXover =  RestXover2(Indd1,:) ;            fprintf(' \n' );
        %------------------

        GetHyperB.findCycleList_maxScafXover([],GetHyperB.SavePremPair,2);
%         GetHyperB.findCycleList_maxScafXover([],GetHyperB.SavePremPair,2);
        
    case 'Min'
        nwhile=1;    
        while nwhile< 5
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
            if Good==1  || size(whichCyl,1)==2   || Goodcc==1 % && sum(SS)==0
                break;
            end
        end
        Good
        whichCyl
end
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
hSurf= GetHyperB.plotScafR_cylindermodel ;

%             ForLegend=surface(nan, nan,'Tag','test','EdgeColor','interp','FaceColor','none','Tag','DefaultScafH');
            % ForLegend=surface(nan, nan, 'Linestyle', 'none', 'Marker', 'none' );
%                         ForLegend=line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');

            hLg= legend(hSurf,'Click me for instructions','Location','northwest' ) ; hLg.String={' \rmScaffold (index \bf\color[rgb]{0.24,0.15,0.66}Start 5'' \color{black}to \color[rgb]{0.78,0.78,0.18}End 3''\rm\color{black})'};
            hLg.Interpreter='tex';        %latex
            hLg.Orientation='horizontal';
%             ForLegend.Marker='.' ; ForLegend.Marker='none';
            hLg.ButtonDownFcn=@(src,evn)LegendBoxing_Scaf( src,evn,ax );
            hLg.Title.String='Click me for instructions' ;
            hLg.Units='normalized'; %hLg.AutoUpdate ='off';
            hLg.Position=[0.0063 0.9528 0.1569 0.0387];

close(fbar) ;

% \color[rgb]{0,0.5,0}

end

