function plotAllCmaes(optDir, outDir, scale)
% Sophy Oliver 2018/2019
% Plot all in an optclim optimisation output directory

if nargin < 2
    outDir = cd;
end
if nargin < 3
    scale = 'Auto';
end

plotParamTrajectCmaes(optDir, outDir, scale)
plotTracerMisfitCmaes(optDir, outDir)
plotParamPairMisfitsCmaes(optDir, outDir, scale)
% plotParamPairProgress(optDir, outDir, scale)
% makeACkpo4Movie(optDir, outDir)
% 
% whereToPlot = 0; % Depth, or 'Atlantic', 'Pacific', 'Indian'
% makeOptMovie(optDir, whereToPlot, outDir)
% 
% plotSpatialMisfits(optDir, generations, outDir)