//LE FICHIER DOIT ETRE IMPERATIVEMENT EN .XLS ET NON .XLSX
function M = importFile(path)
    M=[];
    [fd, SST, Sheetnames, Sheetpos] = xls_open(path);
    [M, TextInd] = xls_read(fd, Sheetpos(1));
    mclose(fd);
    M(:,1)=[];
    M(1,:)=[];
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

//Tracé du nuage en fonction de la matrice
function tracerNuage(M)
    moyenne = mean(M,"r");
    nbIndiv = size(M,"r");
    nbCarac = size(M,"c");
    mat_cov = (M*M) / nbIndiv;
    mat_cor =  Correlation(M);
    echred = M - ones(ni,1)*moyenne;
    ecti = diag((1)./sqrt(diag(mat_cov)));
    echred = echred*ecti;
    
    [D,U]=bdiag(mat_cor); 
    [c,k]=sort(diag(D));
    c=round((c/sum(c,"r"))*100);
    axes = U(:,k);
    indiv = echred*axes;
    
    e= max(abs(indiv));
    xset("window",0); xbasc();
    xset("font",2,3);
    plotframe([-e,-e,e,e],[2,10,2,10],[%f,%f],["Projection sur le plan","1er Axe","2nd Axe"]);
    plot2d(indiv(:,1),indiv(:,2),0,"000");
    cvar = axes([1,2],:);
    cvar = cvar.*.[1,0];
    cvar = cvar';
    plot2d(cvar(:,1),cvar(:,2),5*ones(1,nbCarac),"000");    
    
endfunction



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
    nbLigne = size(Z,"r");
    nbCol = size(BON,"c");
    M = zeros(nbLigne,nbCol);
    for i= 1 : nbLigne
        for j = 1 : nbCol
            M(i,j)= (Z(i,:)*(BON(:,j)));
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

function [tabCaraCentreRed] = calculCarCentreRed (A)
    tabCara=A;
    tabIndi=A';
    n = size(A,"r");
    p = size(A,"c");
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
endfunction


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

function afficherCercleCor(listePoints, numFigure)
    scf(numFigure);
    clf(numFigure);
    theta=0:0.1:2*%pi;
	//Affichage des cercles
    plot(1*cos(theta),1*sin(theta))
	
	
	//Affichage des coordonnees
	taille_listePoints = size(listePoints,"r");
	for i=1:taille_listePoints 
		point = listePoints(i,:);
		plot(point(1,1), point(1,2) ,'+r','markersize',10)
        xstring(point(1,1), point(1,2),string(i));
	end
    plot(0.5*cos(theta),0.5 *sin(theta))
	xgrid
endfunction


function main(pathFileImport, numFigure)
    M = importFile(pathFileImport);
    mprintf('Tableau de données : ');
    disp(M);
    
    n = size(M,"r");
    mprintf("Nombre d individus: ");
    disp(n);
    
    p = size(M,"c");
    mprintf("Nombre de caractères :");
    disp(p);
    
    tabCaracCentreReduit = calculCarCentreRed (M);
    mprintf('Tableau centré réduit :');
    disp(tabCaracCentreReduit);
    
    matCor = Correlation(tabCaracCentreReduit);
    mprintf('Matrice de corrélation :');
    disp(matCor);
    
    [valspropre, base] = calculValeurVecteurBase (matCor);
    
    [qualite,MaxQualiteI,MaxQualiteJ] = qualiteRepresentation(valspropre);
    mprintf("Qualité de la représentation : ")
    disp(qualite);
    
    valsPropreRetenu = [valspropre(MaxQualiteI); valspropre(MaxQualiteJ)];
    basePlan = [base(:,MaxQualiteI), base(:,MaxQualiteJ)];
    mprintf('Valeurs propres retenues: ');
    disp(valsPropreRetenu);
    mprintf('Base du plan:');
    disp(basePlan);
    
    C = composante(basePlan, tabCaracCentreReduit);
    mprintf('Composantes principales:');
    disp(C);
    
    listePointCaractere = coordonneeCaractere(tabCaracCentreReduit,valsPropreRetenu,basePlan, 1, 2);
    afficherCercleCor(listePointCaractere, numFigure);
    
    
endfunction
