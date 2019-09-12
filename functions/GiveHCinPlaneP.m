function [ XY ] = GiveHCinPlaneP( Choice,SelectCase )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if Choice>7
    XY=[];
    return
end


[ HClattice ] = findHClattice( 1 ,[40 40]) ;
HCIndex=cell(1,4)  ;
HCIndex{1}=[ 134  135  148  149  162  163];
HCIndex{2}=[134  135  148  149  162  163  176  177  190  191 ];
% HCIndex{2}=[ 1 8 15 22 ];

HCIndex{3}=[ 134  135  148  149  162  163  176  177  190  191  204  205  218  219 ];
HCIndex{4}=[ 134  135  136  148  149  150  162  163  164  176  177  178  190  191  192  204  205  206  ];
HCIndex{5}=[ 134  148  162  176];
% HCIndex{5}=[ 1 8 15 22 ];
HCIndex{7} =[104 105 117 118 119 120 130 131 132  133  134  135  143  144  145  146  147  148  149  150  157,...
    158  159  160  161  162  163  164  171  172  173  174  175  176  177  178  186  187  188  189  190  191  201  202  203  204  216  217 ];
if Choice==6  % for honeycomb lattic 'X-HC'
[ HClattice ] = findHClattice( 1 ,[40 40]) ;
try 
    load('CustumHCLattice.mat','CustumHCLattice') ;    
    HCIndex{6}=CustumHCLattice{SelectCase} ;
catch
  HCIndex{6}=[37 38 50:53 64 67 78 81 92 95 106:109 121 122 ];  % default custom shaple. Use the end description to setup. Recommend to use helpCADDOM(7) in Command line.
  
end

% HCIndex{6}=[3 4 17:19 31:34 46:49 61:64 76:78 91:14:203 92:14:204 216:218 229:232 242:245 255:258 269:271 283 284 ]; 

end

%-------decide direction :  fixed Apr 24,2019,  prevent this first cylinder
%in custom HC has wrong direction as AGroup
if Choice==6 
ColPlusRowIndex =ceil(HCIndex{6}/HClattice.nrow) + mod(HCIndex{6},HClattice.nrow) ;
CantBeFirst = mod(ColPlusRowIndex,2)==1 ;
end
% for custom crosssection



XY_temp=[HClattice.HCcenter(HCIndex{Choice},1),HClattice.HCcenter(HCIndex{Choice},2)]  ;

[XY ,~]=sortrows(XY_temp);
if Choice==6
    
    XY= circshift(XY ,1) ;
    while 1
        XY= circshift(XY ,-1) ;
        [~, mapToOri] = ismember(XY , XY_temp ,'rows') ;
        if  CantBeFirst(mapToOri(1) )==0
            break
        end
    end
end


%--------example of making more crosssection
%-------   uncommnet the following line. Copy to command line. Write down
%desired cylinder indexes for honeycomb lattice. 
% This example is how to setup default custom HC shape. 
% To change the custom shape, Use helpCADDOM(7) to call GUI to assign. 
% [ HClattice ] = findHClattice( 1 ,[40 40]) ;
% figure(235); scatter(HClattice.HCcenter(:,1) ,HClattice.HCcenter(:,2),'o','filled');
% text(HClattice.HCcenter(:,1)+0.2 ,HClattice.HCcenter(:,2),num2str([1:length(HClattice.HCcenter(:,2))]' ) );
% axis equal ;


%--------

end

