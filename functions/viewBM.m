
    function viewBM(T,BelongTransM)
        if min(BelongTransM)==0
            fprintf('Still some remained unassigned  \n ' )
        end
        figure(4444);clf;hold on;
        for k=1:max(BelongTransM)
        Inds= BelongTransM==k;
        scatter3(T(Inds,1),T(Inds,2),T(Inds,3),'.' )
        end
        axis equal;
    end