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
