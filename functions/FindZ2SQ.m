function [ P ] = FindZ2SQ( Cin,Cout,coor,OriH,CinMoveUp )
%  for SQ
if CinMoveUp==2
    CinMoveUp=0;
end

if OriH<=0   %in case of input smaller than 0
    compentsate=32;
    OriH=OriH+compentsate;
else
    compentsate=0;
end


%%%          (subjectCylinder,nextCylinder,Cylinder coodinate, InputHeight
%%%          by solid model, Up|down,subjectCylindertype)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% XoverPoiInZ=[2 3 4 5 7 8 10 11 12 13 15 16 18 19 20 21 23 24 26 27 28 29 31 32];
% if ~ismember(Cin,ObjAGroup)
if CinMoveUp==0

  G4=[  4 5 15 16 26 27]-ones(1,6) ; 
  G2=[ 10 11 20 21 31 32]-ones(1,6);
  G1=[ 2 3 12 13 23 24]-ones(1,6); 
  G3=[ 7 8 18 19 28 29]-ones(1,6);
elseif CinMoveUp==1
 G2=[  4 5 15 16 26 27];
G4=[ 10 11 20 21 31 32];
 
G3=[ 2 3 12 13 23 24]; 
G1=[ 7 8 18 19 28 29];

   
end


Vec=zeros(1,2);
Vec(1)=int8(coor(Cout,1)-coor(Cin,1));
Vec(2)=int8(coor(Cout,2)-coor(Cin,2));
  Vec=Vec/norm(Vec);
VV=[];
if (Vec(2)==1) &&  (Vec(1)==0)         %North
   VV=1     ;
end
if (Vec(2)==0) &&  (Vec(1)==1)         %East
   VV=2     ;
end
if (Vec(2)==-1) &&  (Vec(1)==0)         %South
   VV=3     ;
end
if (Vec(2)==0) &&  (Vec(1)==-1)         %West
   VV=4     ;
end
Qu=floor(OriH/32);
remain=OriH-Qu*32;
QR=floor(remain/11);
switch CinMoveUp
    case 1                         %up edge 
    switch VV
        case 1             %North
            if   1==1    % ismember(Cin,ObjAGroup)
                 P= Qu*32+G1(QR*2+1);       % P= Qu*32+G1(QR*2-1);    
            else      
                P= Qu*32+G3(QR*2+1);
            end
        case 2             %East
             if  1==1     % ismember(Cin,ObjAGroup)
                 P= Qu*32+G2(QR*2+1);
             else      
                P= Qu*32+G4(QR*2+1);
             end
    
        case 3             %South
             if  1==1     % ismember(Cin,ObjAGroup)
                 P= Qu*32+G3(QR*2+1);
             else      
                P= Qu*32+G1(QR*2+1);
             end
        
        case 4             %West
             if  1==1     % ismember(Cin,ObjAGroup)
                 P= Qu*32+G4(QR*2+1);
             else      
                P= Qu*32+G2(QR*2+1);
             end                        
    end
    
    case 0                        %down
        
    switch VV
        case 1             %North
           if  1==1     %ismember(Cin,ObjAGroup)
                 P= Qu*32+G1(QR*2+2);
            else      
                P= Qu*32+G3(QR*2+2);
            end
            
        case 2             %East
            if  1==1     %ismember(Cin,ObjAGroup)
                 P= Qu*32+G2(QR*2+2);
            else      
                P= Qu*32+G4(QR*2+2);
            end
        case 3             %South
            if 1==1     % ismember(Cin,ObjAGroup)
                 P= Qu*32+G3(QR*2+2);
            else      
                P= Qu*32+G1(QR*2+2);
            end
        
        case 4             %West
           if 1==1     % ismember(Cin,ObjAGroup)
                 P= Qu*32+G4(QR*2+2);
            else      
                P= Qu*32+G2(QR*2+2);
           end              
    end          
end

if compentsate==32
    P=-compentsate+P;
end


end

