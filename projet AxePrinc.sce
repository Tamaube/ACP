
function [](composantes)

    composantes(:,1) = [sqrt(3)/2; sqrt(3)/2; -sqrt(3)]
    composantes(:,2) = [sqrt(3)/2; -sqrt(3)/2; 0]
    disp(composantes)
    
    n=3
    
    lambda = [3/2;1/2]
    
    
    
    for i=1:2
        normeCarComposante = norm(composantes(:,i)) ^2;
        
        if(max(abs(C(i)-sqrt(n*lambda(i))))==0) then
            disp("C ok");
        end
    
    end

endfunction