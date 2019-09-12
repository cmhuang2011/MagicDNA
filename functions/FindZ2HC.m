function [ P ] = FindZ2HC( Cin,Cout,coor,OriH,CinMoveUp )

if CinMoveUp==2
    CinMoveUp=0;
end

if OriH<=0   %in case of input smaller than 0
    compentsate=21;
    OriH=OriH+compentsate;
else
    compentsate=0;
end

Qu=floor(OriH/21);
remain=OriH-Qu*21;
QR=floor(remain/11);

Vec=coor(Cout,:)-coor(Cin,:);
theat=180/pi*atan2(Vec(2),Vec(1));
epson=0.1;

if CinMoveUp==0
  
  if abs(theat-90)<=epson   
   G90=[8  18 ];    
   P= Qu*21+G90(QR+1);       % P= Qu*32+G1(QR*2-1);       
      
  elseif   abs(theat+150)<=epson   %
      G210=[1   11 ];
   P= Qu*21+G210(QR+1);
  elseif   abs(theat+30)<=epson    %
      G330=[4  15 ];
    P= Qu*21+G330(QR+1);   
  end  
  
  
elseif CinMoveUp==1
    
    if    abs(theat-30)<=epson
      G210=[1  11 ];  
       P= Qu*21+G210(QR+1);    
    elseif   abs(theat-150)<=epson    %
       G150=[4  15 ];
       P= Qu*21+G150(QR+1);    
    elseif   abs(theat-270+360)<=epson
       G90=[8  18 ]; 
      P= Qu*21+G90(QR+1);      
    end    
end



% sdfsfsf=34



% if abs(theat-90)<=epson
%     VV=1;
% elseif   abs(theat-210+360)<=epson
%     VV=2;
% elseif   abs(theat-330+360)<=epson
%     VV=3;
% elseif   abs(theat-30)<=epson
%     VV=1;
% elseif   abs(theat-150)<=epson
%     VV=2;
% elseif   abs(theat-270+360)<=epson
%    VV=3; 
% end
% 
% 
% 
% switch CinMoveUp
%     case 0                         %down edge 
%     switch VV
%         case 1             %90
%        
% %                  P= Qu*21+G90(QR+1);       % P= Qu*32+G1(QR*2-1);    
%     
%         case 2             %210
% 
% %                  P= Qu*21+G210(QR+1);
%  
%     
%         case 3             %330
% %                  P= Qu*21+G330(QR+1);
%        
%     end
%     
%     case 1                        %up
%         
%     switch VV
%         case 1             %30
%             P= Qu*21+G210(QR+1);    
%         case 2             %150
%             P= Qu*21+G330(QR+1);    
%         case 3             %270
%             P= Qu*21+G90(QR+1);    
%          
%     end          
% end

if compentsate==21
    P=-compentsate+P;
end
% P
%  [Cin,Cout,OriH]
%  coor
%  [CinMoveUp,DeterCylCinMoveUpect ]
%  ObjAGroup
end

