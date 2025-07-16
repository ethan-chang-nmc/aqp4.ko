% Code by Ethan to analyze confocal Data
% VAR DECLARATION
sheetc = sheetnames("input_excel_sheet.xlsx");
scounter = numel(sheetc);
dconfocal = cell(1,scounter);
nrmldconf = cell(1,scounter);
baseline = cell(1,scounter);
stder =  cell(1,scounter);
decayf = cell(1,scounter);
halff1  = cell(1,scounter);
halff2  = cell(1,scounter);
posh1 = cell(1,scounter);
posh2 = cell(1,scounter);
posd1 = cell(1,scounter);
posd2 = cell(1,scounter);

% for loop to pull data and store into cells
for i = 1:+1:scounter
    dconfocal{1,i} = readmatrix('input_excel_sheet.xlsx','Sheet',i);
end

% normalization to maximum per trace
for i = 1:+1:scounter
    tempmat = dconfocal {i};
    columncount = size(tempmat,2);
    for k = 2:+1:columncount
        tmpmax = max(tempmat(:,k));
        tempmat(:,k) = tempmat(:,k)./ tmpmax;
    end
    nrmldconf{i} = tempmat;
end

% calculating baseline and 1xstandard error value
for i = 1:+1:scounter
    tempmat = nrmldconf {i};
    tempmat1 = baserange{i};
    columncount = size(tempmat,2);
    holding1 = cell(1,columncount);
    holding2 = cell(1,columncount);
    tempval = 0;
    for j = 2:+1:columncount
        tempa = tempmat1{j};
        startbase = min(tempa);
        endbase = max(tempa);
        while isnan (tempmat (endbase,j))
            endbase = endbase - 1;
        end
        totalvals = endbase - startbase;
        A = zeros(1, startbase);
        counter = 1;
        for k = startbase:+1:endbase
            tempval = tempval + tempmat (k,j);
            A(1,counter) = tempmat (k,j);
            counter = counter+1;
        end
        tempval = tempval ./totalvals;
        holding1 {j} = tempval;
        tempval = 0;
        S = (std(A))./(sqrt(length(A)));
        holding2{j} = S;
    end
    baseline {i} = holding1;
    stder{i} = holding2;
end

% decayF
for i = 1:+1:scounter
    tempmat1 = baseline{i};
    tempmat2 = stder{i};
    holding1 = cell(1, columncount);
    columncount = size(tempmat1,2);
    for j = 2:+1:columncount
        tempval = tempmat1{j} + tempmat2{j};
        holding1{j} = tempval;
    end
    decayf{i} = holding1;
end

% half F calculations
for i = 1:+1:scounter
    tempmat = nrmldconf {i};
    tempmat1 = peaks{i};
    columncount = size(tempmat,2);
    rowcount = size(tempmat,1);
    holding1 = cell(1, columncount);
    holding2 = cell(1, columncount);
    tempval1 = 0;
    tempval2 = 0;
    for j = 2:+1:columncount
        tempa = tempmat1 {j};
        firstp = min(tempa); %first peak position
        secondp = max (tempa); %secondpeak position
        firstval = tempmat (firstp, j);
        secondval = tempmat (secondp, j);
        tempval1 = ((firstval + decayf{i}{1,j})/2); %half f 1st
        holding1{j} = tempval1;
        tempval2 = ((secondval + decayf{i}{1,j})/2); %half f 2nd
        holding2{j} = tempval2;
    end
    halff1{i} = holding1;
    halff2{i} = holding2;
end

% finding position values of halfF and converting from pix to micrometers (um)
for i = 1:+1:scounter
    tempmat = peaks {i};
    tempmat1 = nrmldconf{i};
    temph1 = halff1{i};
    temph2 = halff2{i};
    columncount = size(tempmat1,2);
    rowcount = size(tempmat1,1);
    holding1 = cell(1, columncount);
    holding2 = cell(1, columncount);
    for j = 2:+1:columncount
        % skip lectin calculations
        if mod (j-1,3)== 0
            continue
        end
        tempvalh1 = temph1(1,j);
        tempvalh2 = temph2(1,j);
        tempval1 = tempvalh1 {1};
        tempval2 = tempvalh2 {1};
        tempa = tempmat {j};
        firstppos = min(tempa); % first peak
        firstp = tempmat1(firstppos,1);
        secondppos = max (tempa); % secondpeak
        secondp = tempmat1(secondppos,1);
        tempmath1 = zeros(2,rowcount);
        for k = 1:+1:rowcount
            tempmath1 (1,k) = tempmat1 (k,j);
            tempmath1 (2,k) = tempmat1 (k,1);
        end
        % Find position value when it is less than half f and when the
        % position is to the left/right of first/second peak. For first
        % peak, last value in this array should be the halff position - for
        % second peak, first value in this array should be halff position
        [rowh1,colh1] = find (tempmath1(1,:) <= tempval1 & tempmath1(2,:) <= firstp);
        [rowh2,colh2] = find (tempmath1(1,:) <= tempval2 & tempmath1(2,:) >= secondp);
        tempposh1 = tempmat1(colh1(1,size(colh1,2)),1);
        tempposh2 = tempmat1(colh2(1,1),1);
        holding1{j} = tempposh1;
        holding2{j} = tempposh2;
    end
    posh1{i} = holding1;
    posh2{i} = holding2;
end

% finding position values of decayF and converting to
% micrometers (um)
for i = 1:+1:scounter
    tempmat = peaks{i};
    tempmat1 = nrmldconf{i};
    temph1 = decayf{i};
    columncount = size(tempmat1,2);
    rowcount = size(tempmat1,1);
    holding1 = cell(1, columncount);
    holding2 = cell(1, columncount);
    for j = 2:+1:columncount
        % skip lectin calculations
        if mod (j-1,3)== 0
            continue
        end
        tempvalh1 = temph1(1,j);
        tempval1 = tempvalh1 {1};
        tempa = tempmat {j};
        firstppos = min(tempa); % first peak
        firstp = tempmat1(firstppos,1);
        secondppos = max (tempa); % secondpeak
        secondp = tempmat1(secondppos,1);
        tempmath1 = zeros(2,rowcount);
        for k = 1:+1:rowcount
            tempmath1 (1,k) = tempmat1 (k,j);
            tempmath1 (2,k) = tempmat1 (k,1);
        end
        % Find position value when it is less than half f and when the
        % position is to the left/right of first/second peak. For first
        % peak, last value in this array should be the halff position - for
        % second peak, first value in this array should be halff position
        [rowh1,colh1] = find (tempmath1(1,:) <= tempval1 & tempmath1(2,:) <= firstp);
        [rowh2,colh2] = find (tempmath1(1,:) <= tempval1 & tempmath1(2,:) >= secondp);
        tempposh1 = tempmat1(colh1(1,size(colh1,2)),1);
        tempposh2 = tempmat1(colh2(1,1),1);
        holding1{j} = tempposh1;
        holding2{j} = tempposh2;
    end
    posd1{i} = holding1;
    posd2{i} = holding2;
end