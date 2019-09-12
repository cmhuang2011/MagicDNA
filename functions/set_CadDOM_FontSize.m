function set_CadDOM_FontSize(src,evn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% ss=findall(gcf,'-property','FontSize') ;
ss2=findall(gcf,'-property','FontSize','Type','UIControl') ;
set(ss2,'FontSize',src.Value ) ;
% src.Value


ss3=findall(gcf,'-property','FontSize','Type','uitable') ;
set(ss3,'FontSize',src.Value ) ;


end

