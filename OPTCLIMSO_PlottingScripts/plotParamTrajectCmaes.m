function maxGen = plotParamTrajectCmaes(optDir, outDir, scale)
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
[~, paramVals, ~, ~] = i_readFinalJSON(optDir);

% Get Generation Info
genFields = fieldnames(paramVals);
genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
maxGen = max(genNum);
[~,iOrder] = sort(genNum);

% Get Parameter Info
[parNames, parUnits, parRange, numPars] = i_getParamInfo(paramVals);

% Get CMAES Results from Folders (not from JSON file!!)
[~, paramVals, bestPop, bestGen] = i_getCmaesResults(optDir, genFields(iOrder), parNames);

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
    axp(p) = i_plotParamCMAES(axp(p), 1:maxGen, paramVals.(parNames{p}), bestPop, plotNames{p}, parUnits{p}, parRange(p, :), AxisRange);
end

% Annotate Paramaters Values at Best Misfit
for p = 1:numPars
    axp(p) = i_annotateBest(axp(p), paramVals.(parNames{p}), bestPop, bestGen);
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
function hax = i_plotParamCMAES(hax, x, Y, iBest, plotName, parUnit, parRange, AxisRange)

genSize = length(Y(1,:));
popSize = length(Y(:,1));
yBest = NaN(1,genSize);
yBad = NaN(popSize-1, genSize);
for b = 1:genSize
    igood = 1:popSize == iBest(b);
    yBest(b) = Y(igood,b);
    yBad(:,b) = Y(~igood, b);
end
plot(hax, x, yBest, 'k-')
hold(hax, 'on')
plot(hax, x, yBad, 'k.')
if strcmp(AxisRange, 'Auto')
    set(hax, 'ylim', parRange)
end

xlabel(hax, 'Generation', 'fontsize', 16)
ylabel(hax, parUnit, 'fontsize', 16)
title(hax, plotName, 'fontsize', 20)
set(hax, 'fontsize', 16)

% *************************************************************************
function hax = i_annotateBest(hax, Y, bestPop, bestGen)

val2plot = Y(bestPop(bestGen),bestGen);
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