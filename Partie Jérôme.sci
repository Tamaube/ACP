'LE FICHIER DOIT ETRE IMPERATIVEMENT EN .XLS ET NON .XLSX
function M = importFile(path)
    M=[];
    [fd, SST, Sheetnames, Sheetpos] = xls_open(path);
    [M, TextInd] = xls_read(fd, Sheetpos(1));
    mclose(fd);
    M(:,1)=[];
    M(1,:)=[];
endfunction

'Calcul de la qualité de la représentation(en %) par rapport aux valeurs propres choisies.
function quality = qualityOfRepresentation(valsPropresRestantes, valsPropres)
    total = sum(valsPropres);
    totalRest = sum(valsPropresRestantes);    
    quality = totalRest*100/total;
endfunction

'Tracé du nuage en fonction de la matrice
function tracerNuage(M)
    moyenne = mean(M,"r");
    nbIndiv = size(M,"r");
    nbCarac = size(M,"c");
    mat_cov = (M*M) / nbIndiv;
    mat_cor =  Correlation(M);
    echred = M - ones(ni,1)*moyenne;
    ecti = diag((1)./sqrt(diag(mat_cov)));
    echred = echred*ecti;
    
    D,U]=bdiag(mat_cor); 
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
