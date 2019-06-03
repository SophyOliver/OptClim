function plotCompare(optDirs, outDir, isCmaes, names, cols, scale)
% Plot the misit information for all generations of the optimisation
% Example as in plotAllCompare.m
%
% Sophy Oliver Aug-2018 (added CMAES functionality Apr-2019)


%% Set up figs/axes

fig1 = figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
axp1 = axes(fig1);
hold(axp1, 'on')
p1 = gobjects(1,length(optDirs));

fig2 = figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
p2 = gobjects(1,length(optDirs));
        
%% Get optimisation info

for i = 1:length(optDirs)
    
    optDir = optDirs{i};
    isCma = isCmaes(i);
    
    % Read JSON file
    [SumSqCost, paramVals, ~, ~] = i_readFinalJSON(optDir);

    % Get Generation Info
    genFields = fieldnames(SumSqCost);
    genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
    maxGen = max(genNum);
    [~,iOrder] = sort(genNum);

    % Get Parameter Info
    [parNames, ~, parRange, numPars] = i_getParamInfo(paramVals);

    % Get Progression Info
    if ~isCma
        [loopEval, ~, restartEval, softRestarts, ~] = i_readProgTxt2(optDir, maxGen);
    else
        % Get CMAES Results from Folders (not from JSON file!!)
        [SumSqCost, paramVals, bestPop, ~] = i_getCmaesResults(optDir, genFields(iOrder), parNames);
    end

    twinp = [170, 24, 0.03125, 2, 3.2, 0.858];

    %% Misfits

    % Get Cost Values
    costs = NaN(1, maxGen);
    for g = 1:maxGen
        if ~isCma
            ig = genNum == g;
            costs(g) = SumSqCost.(genFields{ig});
        else
            popSize = length(SumSqCost(:,1));
            igood = 1:popSize == bestPop(g);
            costs(g) = SumSqCost(igood,g);
        end
    end

    %% Cost Plot
    
    if ~isCma
        % Reshape for specific plotting requirements
        loopEval(loopEval == 1) = [];
        restartEval = [1, restartEval]; %#ok<AGROW>

        gens = -max(restartEval):max(loopEval)-max(restartEval);
        plot(axp1, gens(1:length(restartEval)), costs(restartEval), [cols{i},'+'], 'linewidth', 2, 'markersize', 8)
        p1(i) = plot(axp1, gens(gens>0), costs(loopEval), [cols{i},'-'], 'linewidth', 2);

        % Plot soft restarts
        i_plotSoft(axp1, softRestarts)
    else
        p1(i) = plot(axp1, 1:maxGen, costs, [cols{i},'-'], 'linewidth', 2);
    end

    %% Params Plot
    
    % Set Up Axes

    if i == 1
        axp2 = gobjects(1,numPars);

        [maxRow, maxCol] = i_getSubplotLayout(numPars);

        for p = 1:numPars
            axp2(p) = subplot(maxRow,maxCol,p);
        end
    end

    % Plot

    plotNames = {'R_O_2_:_P', 'I_C', 'K_P_H_Y', '\muZOO', 'kZOO', 'b*'};
    for p = 1:numPars
        if ~isCma
            [axp2(p), p2(i)] = i_plotParamTraject(axp2(p), parNames{p}, parRange(p, :), paramVals, genFields, genNum, loopEval, restartEval, softRestarts, plotNames{p}, scale, cols{i}, gens);
        else
            thispar = paramVals.(parNames{p});
            par2plot = NaN(1, maxGen);
            for g = 1:maxGen
                igood = 1:popSize == bestPop(g);
                par2plot(g) = thispar(igood,g);
            end
            p2(i) = plot(axp2(p), 1:maxGen, par2plot, [cols{i},'-'], 'linewidth', 2);
        end
        if p > 3
            xlabel(axp2(p), 'Generation')
        else
            set(axp2(p), 'Xtick', [])
        end
    end
    
end

%% Label and Save Figures

% Misfit Plot
legend(axp1, p1, names, 'location', 'northeast', 'fontsize', 30)
xlabel(axp1, 'Generation')
ylabel(axp1, 'Misfit')
set(axp1, 'fontsize', 30)
set(axp1, 'YAxisLocation', 'origin')
%set(axp, 'XLIM', [-13 93])
saveas(fig1, fullfile(outDir, 'plotCompareMisfit.jpg'))

% Param Plot

% Annotate TWIN Parameter Values
for p = 1:numPars
    axp2(p) = i_annotateTWIN(axp2(p), twinp, p);
end

%legend(axp2(2), [p1, p2], names, 'location', 'north', 'orientation', 'horizontal', 'fontsize', 20)

% Annotate Figure Letter
letters = 'abcdef';
for p = 1:numPars
    axp2(p) = i_annotateLetter(axp2(p), letters(p));
end
% set(axp2, 'XLIM', [-13 93])
leg = legend(axp2(p), p2, names, 'orientation', 'horizontal', 'fontsize', 30);
set(leg, 'position', [0.3 0.51 0.4 0.05])
saveas(fig2, fullfile(outDir, 'plotCompareParams.jpg'))

close all
    
%% SUBFUNCTIONS

% *************************************************************************
function [misfits, paramVals, bestPop, bestEval] = i_getCmaesResults(optDir, genFields, parNames)

% Get dimensions info
jname = dir([optDir, '/*final.json']);
jfile = fullfile(optDir, jname(1).name);
fid = fopen(jfile);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
json = jsondecode(str);
popSize = json.optimise.CMAES_namedSettings.populationSize;
genSize = length(genFields);
parSize = length(parNames);

% Set up arrays to fill
misfits = NaN(popSize, genSize);
paramVals = struct;
for p = 1:parSize
    paramVals.(parNames{p}) = NaN(popSize, genSize);
end
bestPop = NaN(1,genSize);

% Misfit info
for n = 1:genSize
	mFile = fullfile(optDir, genFields{n}, 'best_cmaes.txt');
    mVals = csvread(mFile);
    bestPop(n) = mVals(2);
    misfits(:,n) = mVals(3:end);
end

% Params info
for n = 1:genSize
    for d = 1:popSize
        pFile = fullfile(optDir, genFields{n}, ['parameters_input_', num2str(n), '_', num2str(d), '.txt']);
        pVals = csvread(pFile);
        for p = 1:parSize
            paramVals.(parNames{p})(d,n) = pVals(p);
        end
    end
end

% Best info
[~, bestEval] = find(misfits == min(misfits(:)));

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
function i_plotSoft(hax, softRestarts)

yrange = get(hax, 'ylim');
for s = 1:numel(softRestarts)
    plot(hax, [softRestarts(s), softRestarts(s)], yrange, 'r:')
end

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
function [parNames, parUnits, parRange, numPars] = i_getParamInfo(paramVals)

genFields = fieldnames(paramVals);
parNames = fieldnames(paramVals.(genFields{1}));
numPars = numel(parNames);

set6pars = {'ro2ut', 'ACik', 'ACkpo4', 'ACmuzoo', 'AComniz', 'detmartin'};
set6units = {'mmol O_2 mmol^-^1', 'W m^-^2', 'mmol P m^-^3', ...
             '1 d^-^1', '1/(mmol P m^-^3)d^-^1', ''};
set6range = [   150, 200;...
                4.0, 48;...
                0.0001, 0.5;...
                0.1, 4.0;...
                0, 10;...
                0.4, 1.8];
         
[isOK, ipar] = ismember(parNames, set6pars);
if sum(isOK) ~= numPars
    error('Error in plotParamTraject: Code not written to handle units of one/some of these parameters')
end
parUnits = set6units(ipar);
parRange = set6range(ipar, :);

% *************************************************************************
function [maxRow, maxCol] = i_getSubplotLayout(numPars)

if numPars <= 3
    maxRow = 1; maxCol = numPars;
elseif numPars ==4
    maxRow = 2; maxCol = 2;
elseif numPars <=6
    maxRow = 2; maxCol = 3;
else
    error('Error in plotParamTraject: code not written to handle plotting this number of parameters')
end

% *************************************************************************
function [hax, p] = i_plotParamTraject(hax, parName, parRange, paramVals, genFields, genNum, loopEval, restartEval, softRestarts, plotName, AxisRange, Col, gens)

nt1 = length(loopEval);
nt2 = length(restartEval);
vals2plot1 = NaN(1, nt1);
vals2plot2 = NaN(1, nt2);

for g = 1:nt1
    fgood = genNum == loopEval(g);
    thisGen = paramVals.(genFields{fgood});
    parFields = fieldnames(thisGen);
    [~, pgood] = ismember(parName, parFields);
    vals2plot1(g) = thisGen.(parFields{pgood});
end

for g = 1:nt2
    fgood = genNum == restartEval(g);
    thisGen = paramVals.(genFields{fgood});
    parFields = fieldnames(thisGen);
    [~, pgood] = ismember(parName, parFields);
    vals2plot2(g) = thisGen.(parFields{pgood});
end

p = plot(hax, gens(gens>0), vals2plot1, [Col,'-'], 'linewidth', 2);
hold(hax, 'on')
plot(hax, gens(1:length(restartEval)), vals2plot2, [Col,'x'], 'linewidth', 2, 'markersize', 5)

if strcmp(AxisRange, 'Auto')
    set(hax, 'ylim', parRange)
end

yrange = get(hax, 'ylim');
for s = 1:numel(softRestarts)
    plot(hax, [softRestarts(s), softRestarts(s)], yrange, [Col,':'])
end

% y = ylabel(hax, parUnit);
% pos = get(y, 'position');
% rot = get(y, 'rotation');
%set(hax, 'YAxisLocation', 'origin')
% set(y, 'position', pos, 'rotation', rot)
set(hax, 'fontsize', 30)
title(hax, plotName, 'fontsize', 30)
%set(hax, 'xlim', [-15 80])

% *************************************************************************
function hax = i_annotateTWIN(hax, pvals, thisp)

twinVal = pvals(thisp);
hold(hax, 'on')
xrange = get(hax, 'xlim');
plot(hax, xrange, [twinVal, twinVal], 'k--', 'linewidth', 2)

% *************************************************************************
function hax = i_annotateLetter(hax, letter)

pos = get(hax, 'position');
posx = pos(1) - (pos(3)*0.2) + 0.04;
posy = pos(2) + (pos(4)*0.98) - 0.01;
widx = pos(3)*0.05;
widy = pos(4)*0.05;
annotation('textbox', 'String', ['(', letter, ')'],...
    'Position', [posx, posy, widx, widy], 'FontSize', 22, 'LineStyle', 'none');