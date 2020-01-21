function PlotScaf_Fcn(src,evn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ax= findobj(gcf,'Tag','MechScaffold3D');
axes(ax); cltab ; axes(ax);
ss_Assembly= findobj(gcf,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;

%-------------
answer = questdlg('What kind of scaffold routing to plot?', ...
	'Scaffold routing', ...
	'Current','from MagicDNA','from JSON','from JSON');

answer2 = questdlg('How to plot the routing?', ...
	'Scaffold routing', ...
	'Gradient','IsoColor','IsoColor');


% Handle response
switch answer
    case 'Current'
        GetHyperB.plotScafR_cylindermodelMulti(1 ,answer2) ;
%        GetHyperB.plotScafR_cylindermodelMulti(1,) ;

       
    case 'from MagicDNA'
        GetHyperB.plotScafR_cylindermodelMulti(2 ,answer2) ;   
        N= length(GetHyperB.Scaf_fromCadDOM );
    case 'from JSON'
        GetHyperB.plotScafR_cylindermodelMulti(3 ,answer2 ) ; 
        N= length(GetHyperB.Scaf_fromJSON );     
end
N= length(GetHyperB.ScafRouting) ;
str =num2str([1:N]') ;
str=strcat('Scaf-',str) ;
lg= findobj(gctab,'Type','Legend');
lg.String=str ;


end

