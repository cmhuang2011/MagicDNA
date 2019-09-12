% this funcion is used to convert rotation axis representation to rotation
% matrix. In CadDOM, the rotation axis is local axis and convert to global
% rotational matrix and notice always rotate along GC
function A = RotationAxisToRotaionMatrix( s,phi )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
s= s/norm(s) ;
S_M = [ 0 ,-s(3) , s(2) ;  s(3), 0 , -s(1) ; -s(2) , s(1) ,0] ;
A = eye(3) + sind(phi) *S_M + (1- cosd(phi))*S_M*S_M ;


end

