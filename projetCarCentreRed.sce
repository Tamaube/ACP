function [tabCaraCentreRed] = calculCarCentreRed (A)
    tabCara=A;
    tabIndi=A';
    x=0
    //tableau caractère centré réduit 
    tabCaraCentre=zeros(n,p)
    for j=1:p
        for i=1:n
            x = tabCara(i,j) + x
        end
        x = x/n
        for i=1:n
            tabCaraCentre(i,j)=tabCara(i,j)-x
        end  
        x=0
    end
    
    //tableau variance caractère centré
    tabVarCarCentre=zeros(p)
    a = zeros(n)
    b = 0
    for j=1:p
        a = tabCaraCentre(:,j)
        for i=1:n
            b = a(i)^2 + b
        end
        b = b/n
        b = sqrt(b)
        tabVarCaraCentre(j) = b
        b = 0
    end
    //tableau caractère centré Réduit
    tabCaraCentreRed=zeros(n,p)
    
    for j=1:p
    	a = tabCaraCentre(:,j);
        x = a / tabVarCaraCentre(j);
        //x = x * tabVarCaraCentre(j)
        for i=1:n
            tabCaraCentreRed(i,j) = x(i)
        end
    end
    
    //tableau variance caractère centré
    //disp(tabVarCaraCentre)
    //tableau caractère centré réduit 
    //disp(tabCaraCentre)
    //tableau caractère centré réduit 
    disp(tabCaraCentreRed)
endfunction
