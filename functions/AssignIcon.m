function  AssignIcon( btn,imagefilestr )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Locations= strcat(pwd,filesep,'images', filesep) ;
btn.Units='pixels' ;
[a,~]=imread(strcat(Locations,imagefilestr)) ; a=imresize(a, fliplr(btn.Position(1,3:4)));
a(1:2,:,:)=150 ; a(end-1:end,:,:)=150 ; a(:,1:2,:)=150 ; a(:,end-1:end,:)=150 ;
btn.Units='normalized'; btn.CData= a ; btn.String='';
btn.UserData.OriCData=  a ;

% hello=1
end

