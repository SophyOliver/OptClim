function plotTracerMisfit(optDir, outDir)
% Plot the misit information for all generations of the optimisation
%
% optDir = string containing the output directory containing all the
% folders of the optimisation experiment run by OptClimV2
% outDir (optional) = output directory for figure to be saved
%
% Examples: 
% plotTracerMisfit('/Users/orie3677/Documents/Uni/Project/Optimization/OptClim/OxfordMOPS/test_MOPS_bobyqa')
% plotTracerMisfit('/Users/orie3677/Documents/Uni/Project/Optimization/OptClim/OxfordMOPS/test_MOPS_bobyqa', '/Users/outFigures')
%
% Sophy Oliver Aug-2018

%% Get optimisation info

% Read JSON file
[SumSqCost, ~, misfits, bestGen] = i_readFinalJSON(optDir);

% Get Generation Info
genFields = fieldnames(SumSqCost);
genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
maxGen = max(genNum);
bestGen = regexp(bestGen,'\d*','Match');
bestGen = str2double(bestGen{1});

% Get Progression Info
[loopEval, ~, restartEval, softRestarts, ~] = i_readProgTxt2(optDir, maxGen);

%% Misfits

% Get Total Misfit Info
[misNames, tracers, numReg] = i_getMisfitsInfo(misfits.(genFields{1}));

% Get Cost Values
costs = NaN(1, maxGen);
for g = 1:maxGen
    ig = genNum == g;
    costs(g) = SumSqCost.(genFields{ig});
end

% Get Regional Tracer Misfits
RegTrCosts = cell(numel(tracers), 1);

for t = 1:numel(tracers)
    RegTrCosts{t} = NaN(numReg, maxGen);
    thisTr = tracers{t};
    
    for g = 1:maxGen
        ig = genNum == g;
        thisGen = misfits.(genFields{ig});
        
        iTr = cellfun(@(x) strcmp(x(end-2:end), thisTr), misNames);
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

% Calculate Regional Misfits as Fraction of Each Tracer Cost
fracRegMis = cell(numel(tracers), 1);

for t = 1:numel(tracers)
    fracRegMis{t} = RegTrCosts{t}./sumTrMis{t};
end

%% Set Up Axes

fig = figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
axp = gobjects(1,5);

for n = 1:5
    if n <= 2
        axp(n) = subplot(2,3,n);
    else
        axp(n) = subplot(2,3,n+1);
    end
end
pos = get(axp(2), 'Position'); set(axp(2), 'Position', [pos(1) + (pos(3)/1.25), pos(2), pos(3) + (pos(3)/2), pos(4)])
pos = get(axp(1), 'Position'); set(axp(1), 'Position', [pos(1), pos(2), pos(3) + (pos(3)/2), pos(4)])
pos = get(axp(3), 'Position'); set(axp(3), 'Position', [pos(1), pos(2), pos(3), pos(4)])
pos = get(axp(4), 'Position'); set(axp(4), 'Position', [pos(1) - (pos(3)*0.08), pos(2), pos(3), pos(4)])
pos = get(axp(5), 'Position'); set(axp(5), 'Position', [pos(1) - (pos(3)*0.16), pos(2), pos(3), pos(4)])

%% Cost Plot

gens = 1:maxGen;
plot(axp(1), gens(loopEval), log(costs(loopEval)), 'k-')
hold(axp(1), 'on')
plot(axp(1), gens(restartEval), log(costs(restartEval)), 'k.')

% Plot soft restarts
i_plotSoft(axp(1), softRestarts)

% Plot best
i_plotBest(axp(1), bestGen, get(axp(1), 'ylim'))
bestCost = costs(bestGen);
annotation('textbox', 'String', ['Best: ', sprintf('%0.3f', bestCost)],...
    'Position', [0.15, 0.55, 0.1, 0.1], 'FontSize', 22, 'LineStyle', 'none');
xlabel(axp(1), 'Generation')
ylabel(axp(1), 'log (misfit)')

%% Fractional Tracers Plot

tracers = cellfun(@(x) x(~ismember(x, '_')), tracers, 'UniformOutput', false);
nt = numel(tracers);
p = gobjects(1,nt);
c = {'k', 'b', 'r'};

for t = 1:nt
    p(t) = plot(axp(2), gens(loopEval), fracTrMis{t}(loopEval), [c{t}, '-']);
    hold(axp(2), 'on')
    plot(axp(2), gens(restartEval), fracTrMis{t}(restartEval), [c{t}, '.'])
end

% Plot soft restarts
i_plotSoft(axp(2), softRestarts)

% Plot best
i_plotBest(axp(2), bestGen, get(axp(2), 'ylim'))
xlabel(axp(2), 'Generation')
ylabel(axp(2), 'Fraction of Total Misfit')
legend(axp(2), p, tracers, 'location', 'northeast')

%% Fractional Regions Plots

jmap = colormap('jet');
step = floor(length(jmap)/numReg);
nmap = jmap(1:step:numReg*step, :);
upperY = NaN(1, numel(tracers));

for t = 1:numel(tracers)
    for m = 1:numReg
        plot(axp(t+2), gens(loopEval), fracRegMis{t}(m,loopEval),...
            'color', nmap(m,:), 'LineStyle', '-')
        hold(axp(t+2), 'on')
        plot(axp(t+2), gens(restartEval), fracRegMis{t}(m,restartEval),...
            'color', nmap(m,:), 'LineStyle', 'none', 'Marker', '.')
    end
    xlim(axp(t+2), [1 maxGen])
    xlabel(axp(t+2), 'Generation')
    ylabel(axp(t+2), ['Fraction of ', tracers{t}, ' Misfit'])
    colormap(axp(t+2), nmap)
    
    if t == 3
        pos = get(axp(t+2), 'Position');
        c = colorbar(axp(t+2), 'Position', [0.89 pos(2) 0.02 pos(4)]);
        ylabel(c, 'Misfit Region')
        caxis(axp(t+2), [1 20])
    end
    
    yrange = get(axp(t+2), 'Ylim');
    upperY(t) = yrange(2);
    
end

%% Edit Plots

% Re-size axes and plot best and soft restarts
for t = 1:numel(tracers)
    set(axp(t+2), 'Ylim', [0 max(upperY)])
    i_plotBest(axp(t+2), bestGen, get(axp(t+2), 'ylim'))
    i_plotSoft(axp(t+2), softRestarts)
end

% Annotate Figure Letter
letters = 'abcdef';
for p = 1:5
    axp(p) = i_annotateLetter(axp(p), letters(p));
end

%% Save Figure
if nargin ==1
    outDir = cd;
end
saveas(fig, fullfile(outDir, 'plotTracerMisfit.jpg'))

%% Print Information

Header = {'Mean % of Misfit', 'PO4', 'O2', 'NO3', 'Total'};
RowLabels = cell(numReg+1, 1); RowLabels{1} = 'Total';
data = NaN(numReg+1, numel(tracers)+1);

% Total Average Tracer Misfit Fractions
tvals = NaN(1,3);
for t = 1:numel(tracers)
    tvals(t) = mean(fracTrMis{t});
end
data(1,:) = [tvals, sum(tvals)];

% Regional Average Tracer Misfit Fractions
meanFracReg = NaN(numReg, numel(tracers));
for t = 1:numel(tracers)
    meanFracReg(:,t) = mean(fracRegMis{t}, 2);
end

% Regional Average Total Misfit Fractions
fracRegTotMis = (fracRegMis{1}.*fracTrMis{1}) + ...
                (fracRegMis{2}.*fracTrMis{2}) + ...
                (fracRegMis{3}.*fracTrMis{3});
meanFracRegTotMis = mean(fracRegTotMis, 2);   

% Order and place in data
[colTot, isort] = sort(meanFracRegTotMis, 'descend');
colTr = meanFracReg(isort, :);
regs = 1:numReg;
regsort = regs(isort);
for m = 1:numReg
    RowLabels{m+1} = ['Region ', num2str(regsort(m))];
end

data(2:end, 1:end-1) = colTr;
data(2:end, end) = colTot;
data = data.*100;
data = num2cell(data);

disp([Header; [RowLabels, data]])

%% % Plot example data
% rng(1);
% data = NaN(3, maxGen);
% data(1,:) = rand(1, maxGen) * 4 + 7;
% data(2,:) = rand(1, maxGen) * 4 + 3;
% data(3,:) = rand(1, maxGen) * 4 + 2;
% datasum = sum(data);
% datafrac = data ./ datasum;
% figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
% for t = 1:numel(tracers)
%     plot(1:maxGen, datafrac(t, :), 'linewidth', 3)
%     hold on
% end
% xlabel('Generation')
% ylabel('Fraction')
% legend('PO4', 'O2', 'NO3', 'location', 'northeast')
% set(gca, 'fontsize', 28)
% title('Tracer Fraction of Total Misfit', 'fontsize', 42)

% % Plot example data
% figure ('units', 'normalized', 'outerposition', [0 0 1 1]);
% jmap = colormap('jet');
% step = floor(length(jmap)/numReg);
% nmap = jmap(1:step:numReg*step, :);
% t=1;
% data = NaN(numReg, maxGen);
% xmin=0.8;
% xmax=1.2;
% for m = 1:numReg
%     plusby = mean(fracRegMis{t}(m,:)) * 50;
%     if m == 1
%         plusby = plusby*3;
%     end
%     data(m, :) = (xmin+rand(1, maxGen) * (xmax - xmin)) + plusby;
% end
% datasum = sum(data);
% datafrac = data ./ datasum;
% for m = 1:numReg
%     plot(1:maxGen, datafrac(m,:), 'color', nmap(m,:), 'linewidth', 3)
%     hold(gca, 'on')
% end
% xlim(gca, [1 maxGen])
% xlabel(gca, 'Generation')
% ylabel(gca, 'Fraction')
% colormap(gca, nmap)
% c = colorbar;
% ylabel(c, 'Misfit Region')
% set(gca, 'fontsize', 28)
% title('Regional Fraction of PO4 Misfit', 'fontsize', 42)
% caxis(gca, [1 20])
% maxs = max(data,[],2);
% [maxs, imaxs] = sort(maxs);
% regs = [1:19]';
% [regs(imaxs), maxs]

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
function i_plotBest(hax, xVal, yRange)

hold(hax, 'on')
plot(hax, [xVal xVal], yRange, 'k--')


% *************************************************************************
function hax = i_annotateLetter(hax, letter)

pos = get(hax, 'position');
posx = pos(1) - (pos(3)*0.2);
posy = pos(2) + (pos(4)*0.98);
widx = pos(3)*0.05;
widy = pos(4)*0.05;
annotation('textbox', 'String', ['(', letter, ')'],...
    'Position', [posx, posy, widx, widy], 'FontSize', 22, 'LineStyle', 'none');