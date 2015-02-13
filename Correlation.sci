X = fscanfMat("donneesONU.txt");

function COR = Correlation(Z) // fonction qui permet d'avoir la matrice de corrélation des caractères centrés réduits et/ou initiaux
                             //La matrice de corrélation des caractère centrés Yj et Ym est aussi la matrice de covariance des caractères centrés réduits Zj, Zm
	c=size(Z,"c");              //nombre de colonne de la matrice passée en paramètre               
	n=size(Z,"r");              //nombre de ligne de la matrice passée en paramètre
	COR=[]                      // tableau dans lequel on récupère le résultat 
	for i=1:c                   // boucle for qui va de 1 a nbColonne 
		for j=1:c               // idem 
			COR(i,j)=(Z(:,i)'*Z(:,j))/n;    // formule de corrélation avec le produit scalaire afin de récupérer la matrice de corrélation 
		end
	end
endfunction

function M = composante(BON, Z)       // fonction permettant de récuperer toutes les composantes de Z
    M = [];
    nbLigne = size(Z,"r");
    nbCol = size(BON,"c");
    for i= 1 : nbLigne
        for j = 1 : nbCol
            M(i,j)= ((Z(i,:))'*(BON(:,j)));
        end 
    end
endfunction

function Cord = coordonneeCaractere(Z,VPZ,BONZ, composanteI, composanteJ)
    Cord = [];
    nbCol = size(Z,"c");
    for cpt = 1 : nbCol
      Cord = [Cord;[sqrt(VPZ(composanteI)) * BONZ(cpt,composanteI),sqrt(VPZ(composanteJ)) * BONZ(cpt,composanteJ)]]; // on calcule les points (corr(C1,Zm),corr(C2,Zm)) afin de les représenter 
    end
endfunction

function contribution = contributionIndividu(C, lambda) // renvoie la contribution de chaque individu en fonction de chaque composante
    //globales
    nbLigne = size(C,"r");
    nbCol = size(C,"c");
      
    
    for j = 1:nbCol
        somme = 0;
        for i = 1:nbLigne
           somme =somme + C(i,j)*C(i,j); 
        end
        lambda(j)=somme/nbLigne; // calcul des lambda de chaque composant lambda = sigma Cj / n
    end
    contribution = [];
    for i = 1:nbLigne
        somme=0;
        for j= 1:nbCol
           contribution(i,j) = 100*C(i,j)*C(i,j)/(nbLigne*lambda(j)); // contribution en % de l'individu
           somme = somme + contribution(i,j);
        end
        contribution(i,nbCol+1)=somme; // la somme des contributions de tous les individus d'un axe donné devrait être égale à 100%
    end
    for j = 1:nbCol+1
        somme=0;
        for i= 1:nbLigne
           somme = somme + contribution(i,j);
        end
        contribution(nbLigne+1,j)=somme; // la contribution totale d'un individu dans tous les axes
    end
endfunction

//renvoie un vecteur avec les valeurs propres
// et la base du plan principale
function [vals, base] = calculValeurVecteurBase (matCorr)
    [B, D] = spec(matCorr)
    //On trie dans l'ordre décroissant les valeurs propres optenu avec la 
    //fonction spec
    [vals, indices] = gsort(diag(D))
    
    //on recupère le nombre de valeurs dans vals
    taille = length(vals)
    
    // On initialise la variable base
    base = B
    //Comme on a trier les valeurs propres il faut mettre
    //les vecteurs propres associés dans le même ordre
    for  i = 1:taille
        base (:,i) = B(:,indices(i))
    end
endfunction

function afficherCercleCorrelation(listePoints)
	clear
    theta=0:0.1:2*%pi;
	//Affichage des cercles
    plot(1*cos(theta),1*sin(theta))
	plot(0.5*cos(theta),0.5 *sin(theta))
	
	//Affichage des coordonnees
	taille_listePoints = size(listePoints,"c");
	for i:taille_listePoints
		point = listePoints(i)
		plot(point(1), point(2) ,".k")
	end
	xgrid
endfunction

afficherCercles()


mode(0)
clear all

A=[14,10;12,11;10,9]

// probleme 1 
[n,p]=size(A)

tabCara=zeros(n,p)

for i=1:p
    for j=1:n
        tabCara(j,i)=A(j,i)
    end
end

//fin pb1

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


mode(0)
clear all

A=[1/sqrt(2);1/sqrt(2)]
B=[1,1;1,-1]
C=[sqrt(3)/sqrt(2),0,-sqrt(3)/sqrt(2);0,sqrt(3)/sqrt(2),-sqrt(3)/sqrt(2)]


M=zeros(3,2)

for i=1:2
    for j=1:2
        B(i,j)=A(i)*B(i,j)
    end
end

a=0
k=1

for cpt=1:2
    for j=1:3
        for i=1:2
            M(j,cpt) = M(j,cpt) + B(cpt,i) * C(i,k)
        end
        k=k+1
    end
    k=1
end

C=[0;0]

for cpt=1:2
    for j=1:3
        C(cpt)= M(j,cpt)^2 + C(cpt)
    end
    C(cpt) = sqrt(C(cpt))
end
