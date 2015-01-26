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

function afficherCercle(x0, y0, R)
    theta=0:0.1:2*%pi;
    x=x0+R*cos(theta);
    y=y0+R*sin(theta);
    plot(x,y)
endfunction

clear
afficherCercle(0,0,1)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
//afficherCercle(0,0,0.5)
xgrid
