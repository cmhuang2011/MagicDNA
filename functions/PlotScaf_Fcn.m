function PlotScaf_Fcn(src,evn)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
ax= findobj(0,'Tag','MechScaffold3D');
axes(ax); cltab ; axes(ax);
ss_Assembly= findobj(0,'Tag','ss_Assembly') ;
GetHyperB= ss_Assembly.UserData.HyperBundle ;

%-------------
answer = questdlg('What kind of scaffold routing to plot?', ...
	'Scaffold routing', ...
	'Current','from CadDOM','from JSON','from JSON');
% Handle response
switch answer
    case 'Current'
        GetHyperB.plotScafR_cylindermodel(1) ;
    case 'from CadDOM'
        GetHyperB.plotScafR_cylindermodel(2) ;
    case 'from JSON'
        GetHyperB.plotScafR_cylindermodel(3) ;       
end


end

