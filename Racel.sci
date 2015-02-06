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
      Cord = [Cord,[sqrt(VPZ(composanteI)) * BONZ(cpt,composanteI);
	  sqrt(VPZ(composanteJ)) * BONZ(cpt,composanteJ)]]; // on calcule les points (corr(C1,Zm),corr(C2,Zm)) afin de les représenter 
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

