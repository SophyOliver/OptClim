function checkOptMatches(optDir)
% Sophy Oliver 2018/2019
% Check that the model run inputs/outputs match the prog_t6X.txt output
% file and the json file, or highlight any mismatches

%% Read model inputs/outputs

dirOutInfo = i_readDirsInfo(optDir);
    
%% Read text output file

txtOutInfo = i_readProgTxt3(optDir);

%% Read JSON file

jsonOutInfo = i_formatFinalJSON(optDir);

%% Compare

dlen = length(dirOutInfo);
jlen = length(jsonOutInfo);
tlen = length(txtOutInfo);

if dlen == jlen && dlen == tlen
    disp('checkOptMatches: Max eval numbers match')
else
    disp('checkOptMatches: Eval mismatch')
end

% dirVSjson
if dlen == jlen
    dirVSjson = abs(dirOutInfo-jsonOutInfo);
    if max(dirVSjson(:)) < 10^-6
        disp('checkOptMatches: directories and json file data match well')
    else
        disp('checkOptMatches: directiories and json data mismatch')
    end
else
    disp('checkOptMatches: directiories and json eval mismatch')
end

% dirVStxt
if dlen == tlen
    dirVStxt = abs(dirOutInfo-txtOutInfo);
    if max(dirVStxt(:)) < 10^-6
        disp('checkOptMatches: directories and txt file data match well')
    else
        disp('checkOptMatches: directiories and txt data mismatch')
    end
else
    disp('checkOptMatches: directiories and txt eval mismatch')
end

% txtVSjson
if tlen == jlen
    txtVSjson = abs(txtOutInfo-jsonOutInfo);
    if max(txtVSjson(:)) < 10^-6
        disp('checkOptMatches: txt and json file data match well')
    else
        disp('checkOptMatches: txt and json data mismatch')
    end
else
    disp('checkOptMatches: txt and json eval mismatch')
end

% Dir and JSON match but TXT doesn't
if dlen == jlen && tlen ~= dlen
    disp('checkOptMatches: diresctories and json match, but progress text does not')
    
    if tlen < dlen
        disp('checkOptMatches: progress text recorded fewer evaluations')
    else
        disp('checkOptMatches: progress text recorded more evaluations')
        disp('checkOptMatches: Searching for possible duplicates ...')
        [x, ix] = unique(txtOutInfo(:,2));
        if length(x) < tlen
            disp('checkOptMatches: Duplicates found! Duplicated evaluations are:')
            yx = 1:length(txtOutInfo);
            notx = find(~ismember(yx, ix));
            disp(notx)
        else
            disp('checkOptMatches: Duplicates not found')
        end
    end
end
disp('hi')

%% Subfunctions

% *************************************************************************
function [txtOutInfo] = i_readProgTxt3(txtDir)
% Read optimisation progress text

% Read in the progress text
jname = dir([txtDir, '/prog*.txt']);
jfile = fullfile(txtDir, jname(1).name);
fid = fopen(jfile);
tline = fgetl(fid);

txtOutInfo = [];

while ischar(tline)
    if (contains(tline, 'Function eval ') && strfind(tline, 'Function eval ') == 1)
        
        % Eval number
        evalStr = tline(15:20);
        spaces = find(ismember(evalStr, ' '));
        eval = str2double(evalStr(1:spaces(1)));
        
        % F value
        loc1 = strfind(tline, 'f = ') + length('f = ');
        loc2 = strfind(tline, ' at x = ');
        fval = str2double(tline(loc1:loc2-1));
        
        % Parameter values
        loc1 = find(ismember(tline, '['));
        loc2 = find(ismember(tline, ']'));
        
        % What if parameters spread over 2 lines?
        if isempty(loc2)
            pvals1 = tline(loc1+1:end);
            tline = fgetl(fid);
            loc2 = find(ismember(tline, ']'));
            pvals2 = tline(1:loc2-1);
            pvals = [pvals1, pvals2];
        else
            pvals = tline(loc1+1:loc2-1);
        end
        
        % Split by character space
        spaces = ismember(pvals, ' ');
        pvalsm = NaN(1,6);
        for i = 1:6
            idx = find(spaces);
            if i<6
                idx = idx(1) - 1;
                pvalsm(i) = str2double(pvals(1:idx));
                spaces = spaces(idx+1:end); pvals = pvals(idx+1:end);
                idx = find(~spaces); idx = idx(1);
                spaces = spaces(idx:end); pvals = pvals(idx:end);
            else
                if sum(spaces) == 0
                    pvalsm(i) = str2double(pvals);
                else
                    idx = idx(1) - 1;
                    pvalsm(i) = str2double(pvals(1:idx));
                end
            end
        end
        
        % Save
        txtOutInfo(end+1,:) = [eval, fval, pvalsm]; %#ok<AGROW>
        
    end
    tline = fgetl(fid);
end
fclose(fid);

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
function [jsonOutInfo] = i_formatFinalJSON(optDir)
    
% Read JSON file
[costs, paramVals, ~, ~] = i_readFinalJSON(optDir);

% Params
p = {'ro2ut' 'ACik' 'ACkpo4' 'ACmuzoo' 'AComniz' 'detmartin'};

% Generation Info
genFields = fieldnames(paramVals);
genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
[genNum, isort] = sort(genNum);
genFields = genFields(isort);

% jsonOutInfo
eval = genNum;
fval = NaN(size(genNum));
pvals = NaN(length(genNum), length(p));
for i = 1:length(genFields)
    fval(i) = costs.(genFields{i});
    for t = 1:length(p)
        pvals(i,t) = paramVals.(genFields{i}).(p{t});
    end
end
jsonOutInfo = [eval, fval, pvals];

% *************************************************************************
function [dirOutInfo] = i_readDirsInfo(optDir)

dataDir = dir([optDir, '/*dirs']); dataDir = dataDir(1).name;

% Generation Info
genFields = dir(dataDir); genFields=genFields(~ismember({genFields.name},{'.','..'}));
genFields = {genFields.name}';
genNum = cellfun(@(x) str2double(x(end-2:end)), genFields);
[genNum, isort] = sort(genNum);
genFields = genFields(isort);
dirOutInfo = [];

% Loop folders
for i = 1:length(genNum)
    misfits = load(fullfile(dataDir, genFields{i}, 'misfit_output.txt'));
    fval = sum(misfits.^2);
    pvals = load(fullfile(dataDir, genFields{i}, 'parameters_input.txt'))';
    dirOutInfo(end+1,:) = [genNum(i), fval, pvals];
end