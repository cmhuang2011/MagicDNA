function OutP=CheckIntersect(new2Points, oldvertexs,ccase)
if isempty(oldvertexs)
    OutP=[new2Points];
    return
elseif size(oldvertexs,1)<=2
    OutP=[oldvertexs;new2Points];
    return
    
elseif size(oldvertexs,1)>2
    TorF=zeros(size(oldvertexs,1)-2,1);
    for i=1:length(TorF)
       x1=oldvertexs(i,1);y1=oldvertexs(i,2);
       x2=oldvertexs(i+1,1);y2=oldvertexs(i+1,2);
       x3= oldvertexs(end,1);y3= oldvertexs(end,2);
       x4=new2Points(1);y4=new2Points(2);
       
       A=[(y2-y1), -(x2-x1);(y4-y3), -(x4-x3)];
       b=[ (y2-y1)*x1-(x2-x1)*y1 ; (y4-y3)*x3-(x4-x3)*y3 ];
    
       PInt=inv(A)*b;
       if PInt(1)>min(x1,x2) && PInt(1)>min(x3,x4) && PInt(1)<max(x1,x2) && PInt(1)<max(x3,x4)
        TorF(i)=1;
       end             
    end
    if ccase==2;    TorF(1)=0;end
    if sum(TorF)>0  || sum(new2Points==oldvertexs(end,:))==2
        OutP=oldvertexs;
        return
    else
        OutP=[oldvertexs ; new2Points];
        return
    end    
end
end
