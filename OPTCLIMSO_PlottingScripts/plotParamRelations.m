function plotParamRelations(optDirs, outDir, isCmaes, names)
% Plot the misit information for all generations of the optimisation
% Example as in plotAllCompare.m
%
% Sophy Oliver Jan-2019 (added CMAES functionality Apr-2019)

%% Set up figs/axes/outDir/param data

fig = figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
fsize = 14;
axp = gobjects(2,length(optDirs));
p = gobjects(4,length(optDirs));
xranges = NaN(length(optDirs),4);
for i = 1:length(optDirs)
    axp(1,i) = subplot(length(optDirs), 1, i);
    pos = get(axp(1,i), 'position');
    npos = [pos(1)-0.075, pos(2), pos(3)-0.15, pos(4)];
    set(axp(1,i), 'position', npos, 'fontsize', fsize);
    axp(2,i) = axes('Position', npos, 'yaxislocation', 'right', 'color', 'none', 'fontsize', fsize);
    hold(axp(1,i), 'on'); hold(axp(2,i), 'on')
    box(axp(1,i), 'on'); box(axp(2,i), 'on')
    if i == length(optDirs)
        xlabel(axp(1,i), 'Generation');
    end
    title(axp(1,i), names{i})
end

% Set output folder
if ~isdir(fullfile(outDir, 'CompareParameterTuning'))
    mkdir(fullfile(outDir, 'CompareParameterTuning'))
end

% Set parameter sensitivity info (arranged: [PO4, O2, NO3])
parameters = {'ro2ut'; 'ACik'; 'ACkpo4'; 'ACmuzoo'; 'Acomniz'; 'detmartin'};
PO4 = [-0.2661; -0.7056; -0.0028; -0.3935; 0.2094; -0.9823];
O2 = [-0.0183; 2.1214; 0.0252; 0.8784; -0.5624; 4.0486];
NO3 = [0.2844; -1.4157; -0.0224; -0.4849; 0.35; -3.0663];
T = table(PO4, O2, NO3, 'RowNames', parameters);

%% Gen number of figures to make
[~, paramVals, ~, ~] = i_readFinalJSON(optDirs{1});
[parNames, ~, ~, numPars] = i_getParamInfo(paramVals);
[ppairs, ~] = getParamPairs(numPars, parNames);

for f = 1:length(ppairs)
    
    %% Get optimisation info

    for i = 1:length(optDirs)

        optDir = optDirs{i};
        isCma = isCmaes(i);

        % Read JSON file
        [~, paramVals, ~, ~] = i_readFinalJSON(optDir);

        % Get Generation Info
        genFields = fieldnames(paramVals);
        genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
        maxGen = max(genNum);
        [~,iOrder] = sort(genNum);

        % Get Parameter Info
        [parNames, ~, parRange, numPars] = i_getParamInfo(paramVals);

        % Get Parameter Pairs
        [ppairs, ~] = getParamPairs(numPars, parNames);

        % Get Progression Info
        if ~isCma
            [loopEval, ~, restartEval, ~, ~] = i_readProgTxt2(optDir, maxGen);
        else
            % Get CMAES Results from Folders (not from JSON file!!)
            [~, paramVals, bestPop, ~] = i_getCmaesResults(optDir, genFields(iOrder), parNames);
        end

        %% Account for parameter space sampling

        % Reshape for specific plotting requirements
        if ~isCma
            % Reshape for specific plotting requirements
            loopEval(loopEval == 1) = [];
            restartEval = [1, restartEval]; %#ok<AGROW>
            gens = -max(restartEval):max(loopEval)-max(restartEval);
        end

        %% Tracer Plots
        
        if ~isCma
            [vals_Line1, vals_Dash1] = i_ParamValsToPlot(ppairs(f,1), paramVals, genFields, genNum, loopEval, restartEval);
            [vals_Line2, vals_Dash2] = i_ParamValsToPlot(ppairs(f,2), paramVals, genFields, genNum, loopEval, restartEval);
            x_Line = gens(loopEval);
            x_Dash = gens(1:length(restartEval));
        else
            [vals_Line1, vals_Dash1] = i_CmaValsToPlot(paramVals.(ppairs{f,1}), bestPop);
            [vals_Line2, vals_Dash2] = i_CmaValsToPlot(paramVals.(ppairs{f,2}), bestPop);
            x_Line = 1:maxGen;
            x_Dash = 1:maxGen;
        end
       
        p(1,i) = plot(axp(1,i), x_Line, vals_Line1, 'b-', 'linewidth', 2);
        p(3,i) = plot(axp(2,i), x_Line, vals_Line2, 'r-', 'linewidth', 2);
        
        if ~isCma
            p(2,i) = plot(axp(1,i), x_Dash, vals_Dash1, 'b.', 'linewidth', 2, 'markersize', 6);
            p(4,i) = plot(axp(2,i), x_Dash, vals_Dash2, 'r.', 'linewidth', 2, 'markersize', 6);
        else
            for g = 1:length(vals_Dash1(:,1))
                p(2,i) = plot(axp(1,i), x_Dash, vals_Dash1(g,:), 'b.', 'linewidth', 2, 'markersize', 6);
                p(4,i) = plot(axp(2,i), x_Dash, vals_Dash2(g,:), 'r.', 'linewidth', 2, 'markersize', 6);
            end
        end

        % Label
        legend(axp(1,i), [p(1,i), p(3,i)], ppairs{f,1}, ppairs{f,2}, 'location', 'southeast')
        ylabel(axp(1,i), ppairs{f,1});
        ylabel(axp(2,i), ppairs{f,2});
        ip1 = ismember(parNames, ppairs{f,1}); ip2 = ismember(parNames, ppairs{f,2});
        set(axp(1,i), 'YLIM', parRange(ip1,:));
        set(axp(2,i), 'YLIM', parRange(ip2,:));
        xranges(i,1:2) = get(axp(1,i), 'XLIM'); xranges(i,3:4) = get(axp(2,i), 'XLIM');
        
    end
        
    % Same ranges
    minx = min(xranges(:,[1,3])); minx = min(minx);
    maxx = max(xranges(:,[2,4])); maxx = max(maxx);
    for i = 1:length(optDirs)
        set(axp(1,i), 'xlim', [minx maxx]); set(axp(2,i), 'xlim', [minx maxx], 'xticklabel', [])
    end

    % Table
    ip = ismember(parNames, ppairs{f,1}) | ismember(parNames, ppairs{f,2});
    figtab = uitable('Data',T{ip,:}, 'fontsize', 19, 'ColumnName',T.Properties.VariableNames,...
    'RowName',T.Properties.RowNames(ip),'Units', 'Normalized', 'Position',[0.73, 0.45, 0.225, 0.125]);

    annotation('textbox', 'string', 'Table of change in tracer misfits with 10% increase in parameter value',...
                          'position', [0.73, 0.55, 0.225, 0.125], ...
                          'fontsize', 19, 'linestyle', 'none', 'horizontalAlignment', 'center')

    % Save then re-set
    saveas(fig, fullfile('CompareParameterTuning', [ppairs{f,1},'_vs_', ppairs{f,2}, '.png']))
    for i = 1:length(optDirs)
        cla(axp(1,i)); cla(axp(2,i));
        set([axp(1,i), axp(2,i)], 'ylimmode', 'auto', 'xlimmode', 'auto')
    end
    delete(figtab)
    
end

close(fig)


%% SUBFUNCTIONS

function [vals_Best, vals_Rest] = i_CmaValsToPlot(data, bestPop)

genSize = length(data(1,:));
popSize = length(data(:,1));
vals_Best = NaN(1,genSize);
vals_Rest = NaN(popSize-1, genSize);
for g = 1:genSize
    igood = 1:popSize == bestPop(g);
    vals_Best(g) = data(igood,g);
    vals_Rest(:,g) = data(~igood,g);
end
    
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
function [pp, labelp] = getParamPairs(numPars, parNames)

pp = cell(numPars, 2);
ip = 1;
labelp = [];
for p1 = 1:numPars-1
    labelp(end+1) = ip; %#ok<AGROW>
    par1 = parNames(p1);
    par2s = parNames(p1+1:end);
    for p2 = 1:length(par2s)
        pp(ip,1) = par1;
        pp(ip,2) = par2s(p2);
        ip = ip+1;
    end
end

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
function [vals_Loop, vals_Restart] = i_ParamValsToPlot(parName, paramVals, genFields, genNum, loopEval, restartEval)

nt1 = length(loopEval);
nt2 = length(restartEval);
vals_Loop = NaN(1, nt1);
vals_Restart = NaN(1, nt2);

for g = 1:nt1
    fgood = genNum == loopEval(g);
    thisGen = paramVals.(genFields{fgood});
    parFields = fieldnames(thisGen);
    [~, pgood] = ismember(parName, parFields);
    vals_Loop(g) = thisGen.(parFields{pgood});
end

for g = 1:nt2
    fgood = genNum == restartEval(g);
    thisGen = paramVals.(genFields{fgood});
    parFields = fieldnames(thisGen);
    [~, pgood] = ismember(parName, parFields);
    vals_Restart(g) = thisGen.(parFields{pgood});
end