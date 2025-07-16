
    tempmat = nrmldconf {1};% change sheet
    tempmat1 = peaks{1};% change sheet
    columncount = size(tempmat,2);
    rowcount = size(tempmat,1);
    exporthold = zeros(200, 18);
    tempval1 = 0;
    tempval2 = 0;
    fpc1 = 1;
    fpc2 = 2;
    for j = 2:+3:columncount
        tempa = tempmat1 {j};
        firstp = min(tempa); % first peak position
        tfirp1 = firstp;
        secondp = max (tempa); % secondpeak position
        tfirp2 = secondp;
        totval = rowcount - secondp;
        for fp1 = 1:+1:firstp
           exporthold (fp1, fpc1) = tempmat (tfirp1, j); 
           tfirp1 = tfirp1 - 1;
        end
        for fp2 = 1:+1:totval
            exporthold (fp2, fpc2) = tempmat (tfirp2, j); 
            tfirp2 = tfirp2 + 1;
        end
        fpc1 = fpc1 + 2;
        fpc2 = fpc2 + 2;
    end
