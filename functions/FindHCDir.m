function [ Dir ] = FindHCDir(PosCyl,NegCyl,coor)
%edited 08092018
%   Detailed explanation goes here

Vec=coor(NegCyl,:)-coor(PosCyl,:);
theat=180/pi*atan2(Vec(2),Vec(1));

epson=0.1;
  if abs(theat+90)<=epson   %%|| abs(theat-90)<=epson
    Dir=1;
  elseif abs(theat-150)<=epson   %% || abs(theat+30)<=epson
    Dir=2  ;
  elseif abs(theat-30)<=epson    %%|| abs(theat+150)<=epson
    Dir=3  ;      
  end


end

