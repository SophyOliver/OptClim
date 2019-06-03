function plotAllCompare
% Sophy Oliver 2018/2019
% Plot all in an optclim optimisation output directory
% Compares the optimisation progress for 3 separate output directories
% (e.g. each optimised by a different optimisation algorithm)

optDirs = {'/Users/orie3677/Documents/Uni/Project/Chapter1/TWIN/t6b';...
            '/Users/orie3677/Documents/Uni/Project/Chapter1/TWIN/t6d';...
            '/Users/orie3677/Documents/Uni/Project/Optimization/OptClim/OxfordMOPS/outputTest'};
outDir = cd;
names = {'BOBYQA', 'DFO-LS', 'CMA-ES'};
isCmaes = [0,0,1];
cols = {'k', 'b', 'r'};
scale = 'Auto';

plotCompare(optDirs, outDir, isCmaes, names, cols, scale)
plotCompareMisfits(optDirs, outDir, isCmaes, names)
plotParamRelations(optDirs, outDir, isCmaes, names)