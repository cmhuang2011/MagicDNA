function   Output= vec3norm_CM( Vec )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
  Amp = sqrt(sum(Vec.^2 ,2) ) ;  
  R= 1./Amp;
  Output =Vec.*repmat(R,1,3);

end

