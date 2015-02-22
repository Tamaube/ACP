0//renvoie un vecteur avec les valeurs propres
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



