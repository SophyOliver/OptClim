function plotAll(optDir, outDir, scale)
% Sophy Oliver 2018/2019
% Plot all in an optclim optimisation output directory

if nargin < 2
    outDir = cd;
end
if nargin < 3
    scale = 'Auto';
end

checkOptMatches(optDir)
plotParamTraject(optDir, outDir, scale)
plotTracerMisfit(optDir, outDir)
plotParamPairMisfits(optDir, outDir, scale)
% plotParamPairProgress(optDir, outDir, scale)
% makeACkpo4Movie(optDir, outDir)
% 
% whereToPlot = 0; % Depth, or 'Atlantic', 'Pacific', 'Indian'
% makeOptMovie(optDir, whereToPlot, outDir)
% 
% plotSpatialMisfits(optDir, generations, outDir)