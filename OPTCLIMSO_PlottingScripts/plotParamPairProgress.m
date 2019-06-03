function plotParamPairProgress(optDir, outDir, scale)
% Plot the misit information for all parameter pairs of the optimisation as
% frames of a movie, into a separate folder
%
% optDir = string containing the output directory containing all the
% folders of the optimisation experiment run by OptClimV2
% outDir (optional) = output directory for figure directory to be saved
%
% Examples: 
% plotParamPairMisfits('/Users/orie3677/Documents/Uni/Project/Optimization/OptClim/OxfordMOPS/test_MOPS_bobyqa')
% plotParamPairMisfits('/Users/orie3677/Documents/Uni/Project/Optimization/OptClim/OxfordMOPS/test_MOPS_bobyqa', '/Users/outFigures')
%
% Sophy Oliver Aug-2018

%% USER EDITS

AxisRange = scale; % Can be 'Fit' or 'Auto'. 
                    % Auto sets the axes to the parameter max range allowed
                    % Fit tightens the axes to the data

TWIN = 'True';
twinp = [170, 24, 0.03125, 2, 3.2, 0.858];

%% MAIN FUNCTION
                    
% Read JSON file
[SumSqCost, paramVals, ~, bestGen] = i_readFinalJSON(optDir);

% Get Generation Info
genFields = fieldnames(SumSqCost);
genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
maxGen = max(genNum);
bestGen = regexp(bestGen,'\d*','Match');
bestGen = str2double(bestGen{1});

% Get Parameter Info
[parNames, ~, parRange, numPars] = i_getParamInfo(paramVals);

% Get Parameter Pairs
[ppairs, labelp] = getParamPairs(numPars, parNames);

% Get Progression Info
[FuncEval, figIdx, ~] = i_readProgTxt(optDir, maxGen);

% Set up Figure
fig = figure ('units', 'normalized', 'outerposition', [0.225 0 0.55 1]);
maxFig = max(figIdx);

% Get max cost for colorbar
maxMisfit = 0;
for g = 1:length(genFields)
    if SumSqCost.(genFields{g}) > maxMisfit
        maxMisfit = SumSqCost.(genFields{g});
    end
end
maxMisfit = ceil(maxMisfit);

for f = 1:maxFig

    if f == 22
        disp('hi')
    end
    % Identify this function evaluation and previous ones to plot
    ibefore = figIdx < f;
    ithis = figIdx == f;
    gbefore = FuncEval(ibefore);
    gthis = FuncEval(ithis);
    
    % Set Up Axes
    numPlots = sum(1:numPars-1);
    axp = gobjects(1,numPlots);

    % Loop Parameter Pairs and Plot Misfits
    for pp = 1:numPlots

        ypar = ppairs{pp,1};
        xpar = ppairs{pp,2};
        [~,i1] = ismember(ypar, parNames);
        [~,i2] = ismember(xpar, parNames);
        yrangeMax = parRange(i1, :);
        xrangeMax = parRange(i2, :);
        tvaly = twinp(i1);
        tvalx = twinp(i2);
        subp = (i2-1)+((i1-1)*(numPars-1));
        axp(pp) = subplot(numPars-1, numPars-1, subp);
        colormap jet
        caxis([0, maxMisfit])
        hold(axp(pp), 'on')

        % Plot previous misfit(s)
        for ng = 1:length(gbefore)
            g = gbefore(ng);

            % Plot x y c value for this generation of parameter pair
            ig = genNum == g;
            thisGenPar = paramVals.(genFields{ig});

            yval = thisGenPar.(ypar);
            xval = thisGenPar.(xpar);
            cval = SumSqCost.(genFields{ig});
            scatter(axp(pp), xval, yval, [], cval, 'filled')
 

        end
        
        % Plot next misfit(s)
        for ng = 1:length(gthis)
            g = gthis(ng);

            % Plot x y c value for this generation of parameter pair
            ig = genNum == g;
            thisGenPar = paramVals.(genFields{ig});

            yval = thisGenPar.(ypar);
            xval = thisGenPar.(xpar);
            cval = SumSqCost.(genFields{ig});
            scatter(axp(pp), xval, yval, [], cval, 'filled', 'markeredgecolor', 'k')

        end

        % Plot start and TWIN target ref values
        ig = genNum == 1;
        thisGenPar = paramVals.(genFields{ig});
        yval = thisGenPar.(ypar);
        xval = thisGenPar.(xpar);
        plot(axp(pp), xval, yval, 'bp', 'markersize', 14)
        if TWIN
            plot(axp(pp), tvalx, tvaly, 'rp', 'markersize', 14)
        end

        % Bring subplots closer together along x dimension
        axis(axp(pp), 'square')
        pos = get(axp(pp), 'Position');
        xStart = pos(1) - 0.025;
        set(axp(pp), 'Position', [xStart, pos(2:4)], 'TickDir', 'out')

        % Set axis limits according to allowed parameter ranges
        if strcmp(AxisRange, 'Auto')
            set(axp(pp), 'ylim', yrangeMax)
            set(axp(pp), 'xlim', xrangeMax)
        end

        % Label left hand side subplots
        if ismember(pp, labelp)
            xlabel(axp(pp), xpar)
            ylabel(axp(pp), ypar)
        else
            set(axp(pp), 'XTickLabel', [])
            set(axp(pp), 'YTickLabel', [])
        end

        % Colorbar right hand side subplots
        if ismember(pp, labelp-1) || pp == numPlots
            pos = get(axp(pp), 'Position');
            c = colorbar(axp(pp), 'Position', [pos(1) + (pos(3)*1.1), pos(2), 0.01, pos(4)]);
            cticks = get(c, 'Ticks');
            Cticks = arrayfun(@(x) sprintf('%.1f',x), cticks,'un',0);
            set(c, 'TickLabels', Cticks)
            ylabel(c, 'misfit')
        end

    end

    % Title
    if length(gthis) > 1
        titleStr = ['Iterations ', num2str(min(gthis)), ' to ', num2str(max(gthis))];
    else
        titleStr = ['Iteration ', num2str(gthis)];
    end
    annotation('textbox', 'String', titleStr,...
    'Position', [0.1 0.1 0.1 0.1], 'FontSize', 50, 'LineStyle', 'none');
    
    % Save Figure
    if nargin ==1
        outDir = cd;
    end
    if ~isdir(fullfile(outDir, 'plotParamPairProgress'))
        mkdir(fullfile(outDir, 'plotParamPairProgress'))
    end
    saveas(fig, fullfile(outDir, 'plotParamPairProgress', ['pp', num2str(f), '.jpg']))

    clf;
    
end

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
function [FuncEval, figIdx, tlines] = i_readProgTxt(txtDir, maxGen)

jname = dir([txtDir, '/prog*.txt']);
jfile = fullfile(txtDir, jname(1).name);
fid = fopen(jfile);
tline = fgetl(fid);
tlines = {};
while ischar(tline)
    if (contains(tline, 'Function eval ') && strfind(tline, 'Function eval ') == 1) || ...
       (contains(tline, 'Initialising (random directions)') && strfind(tline, 'Initialising (random directions)') == 1) || ...
       (contains(tline, 'Beginning main loop') && strfind(tline, 'Beginning main loop') == 1) || ...
       (contains(tline, 'Restarting from finish point ') && strfind(tline, 'Restarting from finish point ') == 1)
        
        tlines(end+1, 1) = {tline}; %#ok<AGROW>
    end
    tline = fgetl(fid);
end
fclose(fid);

FuncEval = 1:maxGen; FuncEval = FuncEval';
figIdx = NaN(size(FuncEval));
ifig = 1;
ifun = 1;
inloop = 0;

for n = 1:length(tlines)
    tline = tlines{n};
    
    if (contains(tline, 'Function eval ') && strfind(tline, 'Function eval ') == 1)
        figIdx(ifun) = ifig;
        ifun = ifun + 1;
    
    elseif (contains(tline, 'Initialising (random directions)') && strfind(tline, 'Initialising (random directions)') == 1) || ...
            (contains(tline, 'Restarting from finish point ') && strfind(tline, 'Restarting from finish point ') == 1)
        inloop = 1;
        
    elseif (contains(tline, 'Beginning main loop') && strfind(tline, 'Beginning main loop') == 1)
        inloop = 0;
        
    end
    
    if inloop == 0
        ifig = ifig + 1;
    end
    
end