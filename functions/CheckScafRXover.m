function [ Good,whichCyl ] = CheckScafRXover( scafR, tol , RegardlessBundles)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if size(scafR,1) ~= 1 && size(scafR,2) ~= 1  % single scaffold check 
    Good=true;
    whichCyl=-1*ones(200,3); k=1;
    
    SSR=sortrows(scafR) ;
    GCyl=unique(scafR(:,1:2),'rows');
    
    % RegardlessBundles=4:11 ;  % hard
    % RegardlessBundles=[] ;  % hard
    
    inds = ismember(GCyl(:,1) , RegardlessBundles) ;
    GCyl(inds, :)=[] ;   % hard code
    nCyl=size(GCyl,1);
    [aa,bb] =  ismember(SSR(:,1:2) ,GCyl ,'rows') ;    % optimized 09/11/2019
    for cyli=1:nCyl
       ZZ= SSR( bb==cyli ,3) ;
%          ZZ=  SSR(ismember(SSR(:,1:2) ,GCyl(cyli,:)  ,'rows') ,3) ;
        Z2=  sort(ZZ);
        for segj=1:2:length(Z2)-1
            if  Z2(segj+1)- Z2(segj) <tol % &&   Z2(segj+1)- Z2(segj)~=1
                Good=false;
                %             whichCyl= [whichCyl ;GCyl(cyli,:) ];
                whichCyl(k,:) = [GCyl(cyli,:),Z2(segj)]  ; k=k+1;
                whichCyl(k,:) = [GCyl(cyli,:),Z2(segj+1)]  ; k=k+1;
            end
        end
    end    
    whichCyl=setdiff(whichCyl,[-1,-1,-1],'rows') ;    
else   % multi scaffold check 
    Good=true;
    whichCyl=-1*ones(200,3); k=1;
    for k=1 : length(scafR)
        SSR=sortrows(scafR{k}) ;
        GCyl=unique(scafR{k}(:,1:2),'rows');
        
        % RegardlessBundles=4:11 ;  % hard
        % RegardlessBundles=[] ;  % hard        
        inds = ismember(GCyl(:,1) , RegardlessBundles) ;
        GCyl(inds, :)=[] ;   % hard code
        nCyl=size(GCyl,1);
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
        end
    end
    whichCyl=setdiff(whichCyl,[-1,-1,-1],'rows') ;    
end

end

