% code to check if halffs are greater than decayf
checking1 = cell(1, scounter);
checking2 = cell(1,scounter);
storingpos1 = cell(1,20);
storingpos2 = cell(1,20);
counters = 0;
counters2 = 0;
pos1 = 1;
pos2 = 1;
for i=1:+1:scounter
    temps1 = decayf{i};
    temps2 = halff1{i};
    temps3 = halff2{i};
    temps4 = baseline{i};
    columncount = size(tempmat,2);
    holds1 = cell(1,columncount);
    holds2 = cell(1,columncount);
    for j = 2:+1:columncount
        if temps2{j} > temps1{j}
            holds1{j}=1;
        else
            holds1{j}=0;
            if mod(j-1,3)~= 0
                counters = counters + 1;
                shval = i;
                rowval = j;
                storing = [i,j];
                storingpos1 {pos1} = storing;
                pos1 = pos1+1;
            end
        end
        if temps3{j} > temps1{j}
            holds2{j}=1;
        else
            holds2{j}=0;
            if mod(j-1,3)~= 0
                counters = counters + 1;
                shval = i;
                rowval = j;
                storing = [i,j];
                storingpos2 {pos2} = storing;
                pos2 = pos2+1;
            end
        end
        if temps2{j}<temps4{j}
           if mod(j-1,3)~= 0
               counters2 = counters2 + 1;
           end
        end
        if temps3{j}<temps4{j}
           if mod(j-1,3)~= 0
               counters2 = counters2 + 1;
           end
        end
    end
    checking1{i} = holds1;
    checking2{i} = holds2;
end