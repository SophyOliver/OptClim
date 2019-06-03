function maxGen = plotParamTraject(optDir, outDir, scale)
% Optimisation Trajectories of all Parameters
%
% optDir = string containing the output directory containing all the
% folders of the optimisation experiment run by OptClimV2
% outDir (optional) = output directory for figure to be saved
%
% Examples: 
% plotParamTraject('/Users/orie3677/Documents/Uni/Project/Optimization/OptClim/OxfordMOPS/test_MOPS_bobyqa')
% plotParamTraject('/Users/orie3677/Documents/Uni/Project/Optimization/OptClim/OxfordMOPS/test_MOPS_bobyqa', '/Users/outFigures')
%
% Sophy Oliver Aug-2018

%% USER EDITS

twinp = [170, 24, 0.03125, 2, 3.2, 0.858];
AxisRange = scale; % Can be 'Fit' or 'Auto'. 
                    % Auto sets the axes to the parameter max range allowed
                    % Fit tightens the axes to the data

%% Get optimisation info

% Read JSON file
[~, paramVals, ~, bestGen] = i_readFinalJSON(optDir);

% Get Generation Info
genFields = fieldnames(paramVals);
genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
maxGen = max(genNum);

% Get Parameter Info
[parNames, parUnits, parRange, numPars] = i_getParamInfo(paramVals);

% Get Progression Info
[loopEval, ~, restartEval, softRestarts, ~] = i_readProgTxt2(optDir, maxGen);

%% Set Up Axes

fig = figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
axp = gobjects(1,numPars);

[maxRow, maxCol] = i_getSubplotLayout(numPars);
    
for p = 1:numPars
    axp(p) = subplot(maxRow,maxCol,p);
end

%% Plot

plotNames = {'R_O_2_:_P', 'I_C', 'K_P_H_Y', '\muZOO', 'kZOO', 'b*'};
for p = 1:numPars
    axp(p) = i_plotParamTraject(axp(p), parNames{p}, parUnits{p}, parRange(p, :), paramVals, genFields, genNum, loopEval, restartEval, softRestarts, plotNames{p}, AxisRange);
end

% Annotate Paramaters Values at Best Misfit
[~, ibest] = ismember(bestGen, genFields);
bestVals = paramVals.(genFields{ibest});
bestGen = genNum(ibest);
for p = 1:numPars
    axp(p) = i_annotateBest(axp(p), parNames{p}, bestVals, bestGen);
end

% Annotate TWIN Parameter Values
for p = 1:numPars
    axp(p) = i_annotateTWIN(axp(p), twinp, p, max(genNum));
end

% Annotate Figure Letter
letters = 'abcdef';
for p = 1:numPars
    axp(p) = i_annotateLetter(axp(p), letters(p));
end

%% Save Figure

if nargin ==1
    outDir = cd;
end
saveas(fig, fullfile(outDir, 'plotParamTraject.jpg'))

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
        if iEval > maxGen
            return
        end
    
    elseif contains(tline, 'Function eval ') && strfind(tline, 'Function eval ') == 1 && inMainLoop == 0
        evalRestarts(end+1) = iEval; %#ok<AGROW>
        iEval = iEval + 1;
        if iEval > maxGen
            return
        end
    
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
function hax = i_plotParamTraject(hax, parName, parUnit, parRange, paramVals, genFields, genNum, loopEval, restartEval, softRestarts, plotName, AxisRange)

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

plot(hax, loopEval, vals2plot1, 'k-')
hold(hax, 'on')
plot(hax, restartEval, vals2plot2, 'k.')
if strcmp(AxisRange, 'Auto')
    set(hax, 'ylim', parRange)
end

yrange = get(hax, 'ylim');
for s = 1:numel(softRestarts)
    plot(hax, [softRestarts(s), softRestarts(s)], yrange, 'r:')
end

xlabel(hax, 'Generation', 'fontsize', 16)
ylabel(hax, parUnit, 'fontsize', 16)
title(hax, plotName, 'fontsize', 20)
set(hax, 'fontsize', 16)

% *************************************************************************
function hax = i_annotateBest(hax, parName, bestVals, bestGen)

parFields = fieldnames(bestVals);
[~, pgood] = ismember(parName, parFields);
val2plot = bestVals.(parFields{pgood});
hold(hax, 'on')
plot(hax, bestGen, val2plot, 'rd', 'markersize', 14, 'MarkerFaceColor', 'r')

pos = get(hax, 'position');
posx = pos(1) + (pos(3)*0.5);
posy = pos(2) + (pos(4)*0.8);
widx = pos(3)*0.6;
widy = pos(4)*0.15;
annotation('textbox', 'String', ['Best: ', sprintf('%0.3f', val2plot)],...
    'Position', [posx, posy, widx, widy], 'FontSize', 22, 'LineStyle', 'none');

% *************************************************************************
function hax = i_annotateTWIN(hax, pvals, thisp, maxGen)

twinVal = pvals(thisp);
hold(hax, 'on')
plot(hax, [0 maxGen], [twinVal, twinVal], 'k--')

% *************************************************************************
function hax = i_annotateLetter(hax, letter)

pos = get(hax, 'position');
posx = pos(1) - (pos(3)*0.2);
posy = pos(2) + (pos(4)*0.98);
widx = pos(3)*0.05;
widy = pos(4)*0.05;
annotation('textbox', 'String', ['(', letter, ')'],...
    'Position', [posx, posy, widx, widy], 'FontSize', 22, 'LineStyle', 'none');