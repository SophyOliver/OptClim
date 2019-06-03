function plotCompareMisfits(optDirs, outDir, isCmaes, names)
% Plot the misit information for all generations of the optimisation
% Example as in plotAllCompare.m
%
% Sophy Oliver Aug-2018 (added CMAES functionality Apr-2019)

%% Set up figs/axes

fig1 = figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
axp1 = gobjects(1,length(optDirs));
for i = 1:length(optDirs)
    axp1(i) = subplot(length(optDirs), 1, i);
    hold(axp1(i), 'on')
    box(axp1(i), 'on')
    set(axp1(i), 'fontsize', 20);
end

%% Get optimisation info

for i = 1:length(optDirs)
    
    optDir = optDirs{i};
    isCma = isCmaes(i);
    
    % Read JSON file
    [SumSqCost, ~, misfits, ~] = i_readFinalJSON(optDir);

    % Get Generation Info
    genFields = fieldnames(SumSqCost);
    genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
    maxGen = max(genNum);

    % Get Progression Info
    if ~isCma
        [loopEval, ~, restartEval, ~, ~] = i_readProgTxt2(optDir, maxGen);
    end
    
    % Get Total Misfit Info
    [misNames, tracers1, numReg] = i_getMisfitsInfo(misfits.(genFields{1}));
    tracers = {'PO4', 'O2', 'NO3'};


    %% Misfits

    % Get Cost Values
    costs = NaN(1, maxGen);
    for g = 1:maxGen
        ig = genNum == g;
        costs(g) = SumSqCost.(genFields{ig});
    end

    % Get Regional Tracer Misfits
    fracTrMis = calcFracRegMis(tracers, numReg, genNum, maxGen, misfits, genFields, misNames);


    %% Account for parameter space sampling

    if ~isCma
        % Reshape for specific plotting requirements
        loopEval(loopEval == 1) = [];
        restartEval = [1, restartEval]; %#ok<AGROW>
        gens = -max(restartEval):max(loopEval)-max(restartEval);
    end

    %% Tracer Plots

    nt = numel(tracers1);
    c = {'k', 'b', 'r'};
    p1 = gobjects(1,nt);

    for t = 1:nt
        if ~isCma
            p1(t) = plot(axp1(i), gens(loopEval), fracTrMis{t}(loopEval), [c{t}, '-'], 'linewidth', 2);
            plot(axp1(i), gens(1:length(restartEval)), fracTrMis{t}(restartEval), [c{t}, '+'], 'linewidth', 2, 'markersize', 8)
        else
            p1(t) = plot(axp1(i), 1:maxGen, fracTrMis{t}, [c{t}, '-'], 'linewidth', 2);
        end
    end

    ylabel(axp1(i), {'% of Total';[names{i},' Misfit']});
    
end

%% Sort and Save Figure

legend(p1, tracers, 'location', 'northeast', 'fontsize', 20)
xlabel(axp1(i), 'Generation');

ranges = NaN(length(optDirs), 2);
for i = 1:length(optDirs)
    ranges(i,:) = get(axp1(i), 'XLIM');
end
minlim = min(ranges(:,1));
maxlim = max(ranges(:,2));
for i = 1:length(optDirs)
    set(axp1(i), 'XLIM', [minlim, maxlim]);
end

saveas(fig1, fullfile(outDir, 'plotCompareMisTracers.jpg'))

%% SUBFUNCTIONS

% *************************************************************************
function [costd, parameters, simObs, bestEval] = i_readFinalJSON(jsonDir)

jname = dir([jsonDir, '/*final.json']);
jfile = fullfile(jsonDir, jname(1).name);
fid = fopen(jfile);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
json = jsondecode(str);

costd = jsondecode(json.costd);
parameters = jsondecode(json.parameters);
simObs = jsondecode(json.simObs);

cfields = fieldnames(costd);
costvals = NaN(1, length(cfields));
for c = 1:length(cfields)
    costvals(c) = costd.(cfields{c});
end
imin = find(costvals == min(costvals)); imin = imin(1); 
bestEval = cfields{imin};

% *************************************************************************
function [evalInLoop, numRestarts, evalRestarts, softRestarts, tlines] = i_readProgTxt2(txtDir, maxGen)
% Read optimisation progress text

% Read in the progress text
jname = dir([txtDir, '/prog*.txt']);
jfile = fullfile(txtDir, jname(1).name);
fid = fopen(jfile);
tline = fgetl(fid);
tlines = {};
while ischar(tline)
    if (contains(tline, 'Function eval ') && strfind(tline, 'Function eval ') == 1) || ...
       (contains(tline, 'Initialising (random directions)') && strfind(tline, 'Initialising (random directions)') == 1) || ...
       (contains(tline, 'Beginning main loop') && strfind(tline, 'Beginning main loop') == 1) || ...
       (contains(tline, 'Restarting from finish point ') && strfind(tline, 'Restarting from finish point ') == 1) || ...
       (contains(tline, 'Soft restart ') && strfind(tline, 'Soft restart ') == 1)
        
        tlines(end+1, 1) = {tline}; %#ok<AGROW>
    end
    tline = fgetl(fid);
end
fclose(fid);

% Separate the in loop evaluations and the initial/restart evaluations
evalInLoop = [];
numRestarts = [];
evalRestarts = [];
softRestarts = [];

iEval = 1;
inMainLoop = 1;
iRestart = 0;

for n = 1:length(tlines)
    tline = tlines{n};
    
    if contains(tline, 'Function eval ') && strfind(tline, 'Function eval ') == 1 && inMainLoop == 1
        evalInLoop(end+1) = iEval; %#ok<AGROW>
        numRestarts(end+1) = iRestart; %#ok<AGROW>
        iEval = iEval + 1;
    
    elseif contains(tline, 'Function eval ') && strfind(tline, 'Function eval ') == 1 && inMainLoop == 0
        evalRestarts(end+1) = iEval; %#ok<AGROW>
        iEval = iEval + 1;
    
    elseif contains(tline, 'Initialising (random directions)') && strfind(tline, 'Initialising (random directions)') == 1
        inMainLoop = 0;
        
    elseif contains(tline, 'Restarting from finish point ') && strfind(tline, 'Restarting from finish point ') == 1
        inMainLoop = 0;
        iRestart = iRestart + 1;
        
    elseif contains(tline, 'Soft restart ') && strfind(tline, 'Soft restart ') == 1
        inMainLoop = 1;
        iRestart = iRestart + 1;
        softRestarts(end+1) = iEval+1; %#ok<AGROW>
        
    elseif contains(tline, 'Beginning main loop') && strfind(tline, 'Beginning main loop') == 1
        inMainLoop = 1;
        
    end
    
end

evalInLoop(evalInLoop>maxGen) = [];
evalRestarts(evalRestarts>maxGen) = [];

% *************************************************************************
function [misNames, tracers, numReg] = i_getMisfitsInfo(misfits)

misNames = fieldnames(misfits);

tracers = cellfun(@(x) x(end-2:end), misNames, 'UniformOutput', false);
tracers = unique(tracers, 'stable');

misTr1 = cellfun(@(x) strcmp(x(end-2:end), tracers{1}), misNames);
misNames1 = misNames(misTr1);
misNums = cellfun(@(x) x(4:5), misNames1, 'UniformOutput', false);
misNums = cellfun(@(x) str2double(x(~ismember(x, '_'))), misNums);
numReg = max(misNums);

% *************************************************************************
function fracTrMis = calcFracRegMis(tracers, numReg, genNum, maxGen, misfits, genFields, misNames)

RegTrCosts = cell(numel(tracers), 1);

for t = 1:numel(tracers)
    RegTrCosts{t} = NaN(numReg, maxGen);
    thisTr = tracers{t};
    
    for g = 1:maxGen
        ig = genNum == g;
        thisGen = misfits.(genFields{ig});
        
        iTr = cellfun(@(x) strcmp(x(end-(length(thisTr)-1):end), thisTr), misNames);
        thisNames = misNames(iTr);
        
        for m = 1:numReg
            RegTrCosts{t}(m,g) = thisGen.(thisNames{m});
        end
    end
end

% Calculate Tracer Misfits Summed Over Region
sumTrMis = cell(numel(tracers), 1);

for t = 1:numel(tracers)
    sumTrMis{t} = sum(RegTrCosts{t});
end

% Calculate Tracer Misfits as a Fraction of the Total Cost
totTrMis = zeros(1, maxGen);
fracTrMis = cell(numel(tracers), 1);

for t = 1:numel(tracers)
    totTrMis = totTrMis + sumTrMis{t};  
end

for t = 1:numel(tracers)
    fracTrMis{t} = sumTrMis{t}./totTrMis;
end
