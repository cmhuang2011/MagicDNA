
        function checkFcn2(src,~,aH)     
         % used to change axis equal/normal with chechbox
        if src.Value==1
            axis(aH,'normal')
%             axis normal; 
        else
             axis(aH,'equal')
%              axis equal;
        end
        end