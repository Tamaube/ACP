// Projet de Modélisation : ACP 
// Groupe de projet : Steeve VINCENT, Kevin VINCHON, Rodolphe CARGNELLO, Swamynathan CANDASSAMY 
// TP 321

mode(0); // Mode affichage
ieee(1); // Affiche warning pour une exception en virgule flottante
clear all;
normeIncertitude = 0.5;
nbFigure=0;
X = fscanfMat("data.txt","%lg");
disp(X,"Data :")
function Z=MatriceCentreReduite(X)
r=size(X,"r");
c=size(X,"c");
Y=[];
Z=[];
for j=1:c
    column=X(:,j);  
    m = mean(column);
    Y=[Y,column-m*ones(r,1)];// on obtient notre caractère centré Yj
    somme=0;
    for i=1:r
        somme=somme+Y(i,j)*Y(i,j);// somme des Yj
    end
    sigma=somme/r; // correspond a la covariance cov(Xj,Xj) c'est à dire la variance
    sigma = sqrt(sigma); // ce qui nous donne donc l'écart type 
    Z=[Z,(Y(:,j)/sigma)]; // afin d'avoir la mtrice de caractère centré réduit 
end
endfunction

function COR = Correlation(Z) // fonction qui permet d'avoir la matrice de corrélation des caractères centrés réduits et/ou initiaux car la
                             // matrice de corrélation des caractère centrés Yj et Ym c'est aussi la matrice de covariance des caractères centré réduits Zj, Zm
c=size(Z,"c");              // on récupère le nombre de colonne de la matrice passé en paramètre               
n=size(Z,"r");              // on récupère la nombre de ligne de la matrice passé en paramètre
COR=[]                      // déclare le tableau ou on récupère le résultat 
for i=1:c                   // boucle for qui va de 1 a nbColonne 
    for j=1:c               // idem 
        COR(i,j)=(Z(:,i)'*Z(:,j))/n;    // applique la formule de corrélation avec le produit scalaire afin de récupérer la matrice de corrélation 
    end
end
endfunction

function [VP,BON]=valeurPropre(X) // fonction qui va permettre de calculer la valeurs propre et le vecteur propre de la matrice passé en paramètre 
    [BON, VP] = spec(X);            // On récupère les vecteurs propre associé a X et les valeurs propre
    VP=diag(VP);                    //étant donné que les valeurs propre = D dans l'opération précedente, on les replaces en colonne via diag
endfunction

function [Q,MaxQualiteI,MaxQualiteJ] = qualiteRepresentation(VP)         // fonction permettant de récuperer dans une matrice carré toutes les qualités de représentation 
                                                                        //entre chaque caractrèe. Récupère les index des composantes qui donneront la meilleur 
                                                                        //représentation
    n=size(VP,"r");
    Q=[];
    MaxQualiteI = 1;
    MaxQualiteJ = 1;
    maxi=0;
    for i=1:n
        for j=1:n
            Q(i,j)=(VP(i)+VP(j))*100/n;         // qualité de représentation pour toutes les valeurs propres 
            if(Q(i,j)>maxi & i~=j) then
                MaxQualiteI = i;
                MaxQualiteJ = j;
                maxi = Q(i,j);
            end 
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

function Cord = coordonneeCaractere(Z,VPZ,BONZ, composanteI, composanteJ)
    Cord = [];
    nbCol = size(Z,"c");
    for cpt = 1 : nbCol
      Cord = [Cord,[sqrt(VPZ(composanteI)) * BONZ(cpt,composanteI);sqrt(VPZ(composanteJ)) * BONZ(cpt,composanteJ)]]; // on calcul les points (corr(C1,Zm),corr(C2,Zm)) afin de le représenter 
    end
endfunction

function Q2 = QualiteRepresentationIndividu(BON,Z, composanteI, composanteJ)     
    // fonction permettant de trouver les points pour faire la projection des individus
    // on cherche Pi = (Zi.u1/norme(Zi);Zi.U2/norme(Zi))
    nbIndividu = size(Z,"r");
    nbAxe = size(BON,"c");
    Q2 = zeros(nbIndividu,2);
    for i = 1 : nbIndividu
          scal = ((Z(i,:)')'*(BON(:,composanteI)));
          norme = norm(Z(i,:));
          Q2(i,1) = (scal*scal/(norme*norme));
          scal = ((Z(i,:)')'*(BON(:,composanteJ)));
          norme = norm(Z(i,:));
          Q2(i,2) = (scal*scal/(norme*norme));
    end
endfunction

function CordProj = coordonneeIndividu(BON,Z, composanteI, composanteJ) // fonction qui calcul les coordonnées pour le nuage de points 
    nbIndividu = size(Z,"r");
    nbAxe = size(BON,"c");
    CordProj = zeros(2,nbIndividu);
    for i = 1 : nbIndividu
          scal = ((Z(i,:)')'*(BON(:,composanteI)));
          norme = norm(Z(i,:));
          CordProj(1,i) = (scal/(norme));
          scal = ((Z(i,:)')'*(BON(:,composanteJ)));
          norme = norm(Z(i,:));
          CordProj(2,i) = (scal/(norme));
    end
endfunction;


function contribution = contributionIndividu(C) // renvoie la contribution de chaque individu en fonction de chaque composantes
    //globales
    nbLigne = size(C,"r");
    nbCol = size(C,"c");
    lambda = [];
    
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
        contribution(i,nbCol+1)=somme; // la somme des contributions de tous les individus d'un axes donné devrait égale à 100%
    end
    for j = 1:nbCol+1
        somme=0;
        for i= 1:nbLigne
           somme = somme + contribution(i,j);
        end
        contribution(nbLigne+1,j)=somme; // la contribution totale  d'un individu dans tous les axes
    end
endfunction


function nuage(Coord, Axei , Axej, numFigure) // Fonction permettant d'avoir le nuage de point 
    scf(numFigure);
    clf(numFigure);
  xset("font",4,3);

  titre="Projection des individus sur le plan  ";

  n = size(Coord,"r");  
  for l=1:n
        //projection des individus dans le plan i-j
    plot(Coord(l,Axei),Coord(l,Axej),'+r','markersize',10);
    xstring(Coord(l,Axei),Coord(l,Axej),string(l));
  end;
  xtitle(titre);
  xgrid;
endfunction;

function CercleCorrelation(CORD, titre, numFigure) // fonction qui permet de tracer le cercle de corrélation et des projections des caractères
    scf(numFigure);
    clf(numFigure);
    R=1;    // on défini le rayon pour le cercle unitaire
    r=0.5;  // in défini le rayon pour le cercle d'incertitude
    x0=0;
    y0=0;
    theta=0:0.2:%pi*%pi; // theta permet de tracer le cercle unitaire 
    x=x0+R*cos(theta);  // on défini l'équation permettant de tracer le cercle
    y=y0+R*sin(theta);
    plot(x,y);          // plot permet tracer concretement le cercle 
    
    
    CordColonne = size(CORD,"c"); // on récupère le nombre de colonne de caractère 
    for i = 1 : CordColonne
        plot(CORD(1,i),CORD(2,i), '+r','markersize',10); // on place les points de coordonnées de caractère 
        xstring(CORD(1,i),CORD(2,i),string(i)); // placer les légendes de chaque points
    end
    X = x0+r*cos(theta);
    Y = y0+r*sin(theta);
    plot(X,Y); // on trace le 2eme cercle (incertitutde)
    
    xtitle(titre);
    
    xgrid;
endfunction

function NouveauX = ACP(DATA, numFigure)
    nbCaractere=size(DATA,"c");
    nbIndividu=size(DATA,"r");
    disp(nbIndividu,"nombre d''individu :")
    disp(nbCaractere,"nombre de caractères :")
    Z = MatriceCentreReduite(DATA);
    disp(Z,"Matrice Centre Reduite=");
    CORX = Correlation(DATA);
    CORZ = Correlation(Z);
    disp(CORZ,"Matrice Correlation des caractères centrés réduits (CORZ)");
    disp(CORX,"Matrice Correlation des caractères initiaux (CORX)");
    [VPX,BONX] = valeurPropre(CORX);
    [VPZ,BONZ] = valeurPropre(CORZ);

    disp(VPX,"Valeurs prorpres de la matrice de corrélation COORX");
    disp(BONX,"Vecteurs propres de la matrice de corrélation COORX");

    disp(VPZ,"Valeurs prorpres de la matrice de corrélation COORZ");
    disp(BONZ,"Vecteurs propres de la matrice de corrélation COORZ");

    [Q,meilleurI,meilleurJ] = qualiteRepresentation(VPZ);
    disp(Q,"Qualité de représentation pour chaque valeur propres de VPZ")
    disp(BONZ(:,meilleurI),"Meilleur axe 1")
    disp(BONZ(:,meilleurJ),"Meilleur axe 2")
    
    C = composante(BONX,Z);
    disp(C,"nouvelles coordonées de chaques caractères (C)")
    disp(C(:,meilleurI),"Coordonnée C1") // Ici on récupère les 2 composantes principales grâce aux traitements effectués précedemment 
    disp(C(:,meilleurJ),"Coordonnée C2")

    cordCaractere = coordonneeCaractere(Z,VPZ,BONZ,meilleurI,meilleurJ); //
    disp(cordCaractere,"Coordonnées des caractères dans le cercle de corrélation (Représentation 2D)")
    QIndividu = QualiteRepresentationIndividu(BONZ,Z,meilleurI,meilleurJ);
    disp(QIndividu,"Qualités de représentation des Individus");
    cordIndividu = coordonneeIndividu(BONZ,Z,meilleurI,meilleurJ);
    
    disp(cordIndividu,"Coordonnée des qualités de projection des individus")
    
    contribution = contributionIndividu(C);
    disp(contribution, "contribution des individus dans chaque composantes")
    nbFigure=3*(numFigure-1)+1;
    CercleCorrelation(cordCaractere, "Projection des caractéres", nbFigure);
    nbFigure=nbFigure+1;
    nuage(C,meilleurI,meilleurJ, nbFigure);
    nbFigure=nbFigure+1;
    CercleCorrelation(cordIndividu, "Projection des Individus", nbFigure);
    Xintermediaire = [];
    for i = 1 : nbIndividu
       norme = norm(cordIndividu(:,i));
       if(norme>normeIncertitude) then // on ne garde que les individus qui sont aux normes
           Xintermediaire = [Xintermediaire;DATA(i,:)]
       end
    end
    NouveauX = [];
    for j = 1 : nbCaractere
       norme = norm(cordCaractere(:,j));
       if(norme>normeIncertitude) then // on ne garde que les caractères qui sont aux normes
           NouveauX = [NouveauX,Xintermediaire(:,j)]
       end
    end
endfunction

NewDATA = ACP(X,1);
disp(NewDATA,"Nouvelles données")
ACP(NewDATA,2); // on réeffectue l'ACP sans les individus et les caractères mal représenté
