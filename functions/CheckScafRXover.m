function [ Good,whichCyl ] = CheckScafRXover( scafR, tol , RegardlessBundles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

SSR=sortrows(scafR) ;
GCyl=unique(scafR(:,1:2),'rows');

% RegardlessBundles=4:11 ;  % hard 
% RegardlessBundles=[] ;  % hard 

inds = ismember(GCyl(:,1) , RegardlessBundles) ;
GCyl(inds, :)=[] ;   % hard code 

nCyl=size(GCyl,1);
Good=true;
whichCyl=-1*ones(200,3); k=1;
for cyli=1:nCyl
    ZZ=  SSR(ismember(SSR(:,1:2),GCyl(cyli,:)  ,'rows'),3) ;
    Z2=sort(ZZ);
    for segj=1:2:length(Z2)-1
        if  Z2(segj+1)- Z2(segj) <tol % &&   Z2(segj+1)- Z2(segj)~=1
            Good=false;
            %             whichCyl= [whichCyl ;GCyl(cyli,:) ];
            whichCyl(k,:) = [GCyl(cyli,:),Z2(segj)]  ; k=k+1;
            whichCyl(k,:) = [GCyl(cyli,:),Z2(segj+1)]  ; k=k+1;
 
        end
    end
%     for segj=1:length(Z2)-1
%         if  Z2(segj+1)- Z2(segj) <tol  &&   Z2(segj+1)- Z2(segj)~=1
%             Good=false;
%             %             whichCyl= [whichCyl ;GCyl(cyli,:) ];
%             whichCyl(k,:) = GCyl(cyli,:) ; k=k+1;
%         end
%     end
end
whichCyl=setdiff(whichCyl,[-1,-1,-1],'rows') ;

% 


end

