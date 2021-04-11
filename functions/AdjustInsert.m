function  AdjustInsert( backboneStrand , centerH )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% input: a PAIR of lines: one is backbone and the other is center. Use in
% oxDNAInitial.m to seperate insertion bases. 

if length(backboneStrand.XData)~=length(centerH.XData)
   fprintf('Warning. Number of bases not matched in center and backbone in AdjustInsert.m !!  \n')
    return
end

XYZ_back= [backboneStrand.XData ; backboneStrand.YData ;backboneStrand.ZData ]' ;
XYZ_Center= [centerH.XData ; centerH.YData ;centerH.ZData ]' ;

% [C,IA,IC] = uniquetol(XYZ_back ,'ByRows',true,'OutputAllIndices',true)

% [aa,bb] = sort(IC)
XYZ_back(:,end+1 ) =0 ; tol=1e-6 ; cc=1 ;

XYZ_back(1,4)=cc ;
for k=1: size(XYZ_back,1)-1
    if norm( XYZ_back(k,1:3) -XYZ_back(k+1,1:3))< tol
       XYZ_back(k+1,4) = cc ;
    else
       XYZ_back(k+1,4) = cc+1 ;
       cc=cc+1;
    end  
end

RepeatArr = XYZ_back(:,4) ;
Diff_Rp = diff(RepeatArr)==0 ;
Dup = [Diff_Rp; 0]+[0 ;Diff_Rp]  ;
WhereToSep = bwlabel(Dup) ;
for k=1:max(WhereToSep)
  BaseInd= find(WhereToSep==k) ;
  if BaseInd(1)==1
      Pinch_Aind=  1 ;
  else
      Pinch_Aind=  BaseInd(1)-1 ; 
  end
  
  if BaseInd(end)==length(WhereToSep)
      Pinch_Bind=  length(WhereToSep) ;
  else
      Pinch_Bind=  BaseInd(end)+1 ; 
  end
  
  OneD_interp = interp1([Pinch_Aind , Pinch_Bind]',[ XYZ_back(Pinch_Aind,1:3) ;XYZ_back(Pinch_Bind,1:3) ],BaseInd) ;
  XYZ_back(BaseInd,1:3 ) = OneD_interp ;
  
  OneD_interp2 = interp1([Pinch_Aind , Pinch_Bind]',[ XYZ_Center(Pinch_Aind,1:3) ;XYZ_Center(Pinch_Bind,1:3) ],BaseInd) ;
  XYZ_Center(BaseInd,1:3 ) = OneD_interp2 ;
 
  
end


backboneStrand.XData= XYZ_back(:,1)' ;
backboneStrand.YData= XYZ_back(:,2)' ;
backboneStrand.ZData= XYZ_back(:,3)' ;

centerH.XData= XYZ_Center(:,1)' ;
centerH.YData= XYZ_Center(:,2)' ;
centerH.ZData= XYZ_Center(:,3)' ;


end

