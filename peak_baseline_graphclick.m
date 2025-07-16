% Code by Ethan for Peak and Baseline Calculator
% VAR DECLARATION
sheetc = sheetnames("input_excel_sheet.xlsx");
scounter = numel(sheetc); % how many sheets
dconfocal = cell(1,scounter); % cell with each sheet
nrmldconf = cell(1,scounter); % cell with each normalized sheet
peaks = cell(1,scounter);
baserange = cell (1,scounter);

% for loop to pull data and store into cells
for i = 1:+1:scounter
    dconfocal{1,i} = readmatrix('input_excel_sheet.xlsx','Sheet',i);
end

% Normalization to maximum per trace
for i = 1:+1:scounter
    tempmat = dconfocal {i};
    columncount = size(tempmat,2);
    for k = 2:+1:columncount
        tmpmax = max(tempmat(:,k));
        tempmat(:,k) = tempmat(:,k)./ tmpmax;
    end
    nrmldconf{i} = tempmat;
end

% Peak X position identification
%% USAGE: CHOOSE PEAKS FIRST, THEN BASELINE RANGE. ALWAYS CLICK LEFT TO RIGHT! %%
for i = 1:+1:scounter
    tempmat1 = nrmldconf{i};
    columncount = size(tempmat1,2);
    holdPeak = cell(1, columncount);
    holdBase = cell (1,columncount);
    for j = 2:+1:columncount
        tempmat2 = tempmat1(:,j);
        % opens graph for peak X identification
        figure
        plot(tempmat2)
        [xp,yp] = ginput(2);
        holdPeak{j} = round(xp);
        close
        % opens graph for baseline X range identification
        figure
        plot(tempmat2)
        [xb,yb] = ginput(2);
        holdBase{j} = round(xb);
        close
    end
    peaks{i} = holdPeak;
    baserange{i} = holdBase;
end