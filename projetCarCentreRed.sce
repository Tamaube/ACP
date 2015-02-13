mode(0)
clear all
function tabCara = Car(X)
//A=[14,10;12,11;10,9]

// probleme 1 
[n,p]=size(X)

tabCara=zeros(n,p)

for i=1:p
    for j=1:n
        tabCara(j,i)=X(j,i)
    end
end
endfunction

//fin pb1

function tabIndi = Ind(X)
//transpose il y a une fonction pour sa
tabIndi=zeros(p,n)

for i=1:p
    for j=1:n
        tabIndi(i,j)=A(j,i)
    end
end

x=0

//good boy !!! mais choisi des nom EXPLICITE
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

	//tu peux eviter a boucle genre a = tabCaraCentre(:,j)
    for i=1:n
        a(i) = tabCaraCentre(i,j)
    end
	
	// la variance c'est la norme au carré divisé par n utilise les fonction de scilab!!!!!!
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
	//tu peu eviter a boucle genre a = tabCaraCentre(:,j)
    for i=1:n
        a(i) = tabCaraCentre(i,j)
    end
	//nom c'est la racine carrée de la variance pas la variance tout cour
    x = a / tabVarCaraCentre(j)
    //x = x * tabVarCaraCentre(j)
    for i=1:n
        tabCaraCentreRed(i,j) = x(i)
    end
end

//tableau variance caractère centré
tabVarCaraCentre
//tableau caractère centré réduit 
tabCaraCentre
//tableau caractère centré réduit 
tabCaraCentreRed

