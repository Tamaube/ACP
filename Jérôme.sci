function M = importFile(path)
    M=[];
    [fd, SST, Sheetnames, Sheetpos] = xls_open(path);
    [M, TextInd] = xls_read(fd, Sheetpos(1));
    mclose(fd);
    M(:,1)=[];
    M(1,:)=[];
endfunction