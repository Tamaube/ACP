function [] = AxePrinc (composantes,lambda)

    for i=1:2
        normeCarComposante = norm(composantes(:,i)) ^2;
        
        if(max(abs(C(i)-sqrt(n*lambda(i))))==0) then
            disp("C ok");
        end
    
    end

endfunction

function M = composante(BON, Z)       // fonction permettant de récuperer toutes les composantes de Z
    M = [];
    nbLigne = size(Z,"r");
    nbCol = size(BON,"c");
    for i= 1 : nbLigne
        for j = 1 : nbCol
            M(i,j)= ((Z(i,:)')'*(BON(:,j)));
        end 
    end
endfunction
