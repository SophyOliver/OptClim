%% Make biomes mask files for MOPS model to calculate misfit for each region

% Set toplevel path to GCMs configuration

base_path='/Users/orie3677/Documents/Original_TMM_MOPS/TransportMatrixConfigs/MITgcm_2.8deg';
%base_path='/Users/orie3677/Documents/Original_TMM_MOPS/TransportMatrixConfigs/MITgcm_ECCO';
%base_path='/Users/orie3677/Documents/Original_TMM_MOPS/TransportMatrixConfigs/UVicKielIncrIsopycDiff';

biomes_file = '/Users/orie3677/Documents/Uni/Project/Chapter1/Biomes/BiomesData.mat';

% Set path names, etc.
writeFiles=1;
load(fullfile(base_path,'config_data'))

matrixPath=fullfile(base_path,matrixPath);

gridFile=fullfile(base_path,'grid');
boxFile=fullfile(matrixPath,'Data','boxes');
profilesFile=fullfile(matrixPath,'Data','profile_data');

load(gridFile,'nx','ny','nz','dznom','dz','da','x','y','z','deltaT','gridType')
load(boxFile,'Xboxnom','Yboxnom','Zboxnom','izBox','nb','volb')

%% Regions for Misfit avarages (Sophy Oliver 12-Jul-18)

% Load biomes
biomes = load(biomes_file);
fnames = fieldnames(biomes);
biomes = biomes.(fnames{1});

xbio = biomes.x;
ybio = biomes.y;
zbio = biomes.z;
biomes = biomes.biomes;

% Re-grid biomes

% Pad to extend the edges of the grid so that interp3 nearest can work at the real edges
xbio = [-9999;xbio;9999];
ybio = [-9999;ybio;9999];
zbio = [-9999;zbio;9999];
xpad = [-9999;x;9999];
ypad = [-9999;y;9999];
zpad = [-9999;z;9999];

[XBIO,YBIO,ZBIO] = meshgrid(xbio,ybio,zbio);
[XGRID,YGRID,ZGRID] = meshgrid(xpad,ypad,zpad);

biomes = padarray(biomes,[1 1 1],0,'both');

% Interp
biomes_ngrid = interp3(XBIO,YBIO,ZBIO,biomes, XGRID,YGRID,ZGRID, 'nearest');
biomes = round(biomes_ngrid);

% Cut off edges
biomes = biomes(2:end-1,2:end-1,2:end-1);
XGRID = XGRID(2:end-1,2:end-1,2:end-1);
YGRID = YGRID(2:end-1,2:end-1,2:end-1);
ZGRID = ZGRID(2:end-1,2:end-1,2:end-1);

% % % Plot biomes to check
% for d = 1:length(z)
%     
%     surf(XGRID(:,:,1),YGRID(:,:,1),zeros(size(XGRID(:,:,1))), biomes(:,:,d), 'edgecolor', 'none')
%     view(0,90)
%     pause
% end

% Create masks for each region

regions = unique(biomes); regions(isnan(regions)) = [];
numRegions = length(regions);

% Create masks
amask = zeros(nb,numRegions);

% Loop boxnom values
for i = 1:nb

    ix = find(x == Xboxnom(i));
    iy = find(y == Yboxnom(i));
    iz = find(z == Zboxnom(i));
    thisreg = biomes(iy,ix, iz);

    if ~isnan(thisreg)
        ireg = find(regions == thisreg);
        amask(i,ireg) = 1;
    end

end

% At the moment we're just duplicating this mask for all 3 observational tracers
PO4amask = amask;
NO3amask = amask;
O2amask = amask;

% % Plot regions to check

% % Plot
% f1 = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
% axesm('robinson', 'maplonlimit', [0 360], 'frame', 'on');
% set(gcf,'renderer','zbuffer')
% contourm(YGRID(:,:,1), XGRID(:,:,1), biomes(:,:,1), 'LineColor', 'k', 'LineWidth', 3, ...
%         'LevelStep', 1.001)
% hold on
% geoshow('landareas.shp', 'facecolor', [0.5 0.5 0.5]);
% for reg = 1:numRegions
%     igood = PO4amask(:,reg) == 1;
%     splot = scatterm(Yboxnom(igood),Xboxnom(igood),25,PO4amask(igood,reg),'filled');
%     minz = min(Zboxnom(igood));
%     maxz = max(Zboxnom(igood));
%     title(['Region ', num2str(reg), '. Depth Range = ', num2str(minz), ' - ', num2str(maxz), 'm'])
%     pause
%     delete(splot)
% end

%% Write files

if writeFiles
  
  % Regional weighting files  
  writePetscBin('volb.petsc',volb)
  writePetscBin('PO4_amask.petsc',PO4amask,1)
  writePetscBin('NO3_amask.petsc',NO3amask,1)
  writePetscBin('O2_amask.petsc',O2amask,1)
  
end