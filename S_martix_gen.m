% Convert the ModelSEED ID to stoichiometric matrix (S)
% Build the model structure .mat file
% ref: https://github.com/ModelSEED/ModelSEEDDatabase
%% load data (revised based on the file name)
% Import the Excel data into MATLAB
% Assuming 'GSM.xlsx' are the Excel files.
reactions = readtable('GSM.xlsx',sheet = 'iCAC802 Reactions');
metabolites = readtable('GSM.xlsx',sheet = 'iCAC802 Metabolites');

%% generate S matrix
% Initialize a map to keep track of metabolite indices
reactions = reactions.SEED_FORMULA;  % e.g., cpd00035 + cpd00040 <==> cpd00020 + cpd00033
metaboliteList =metabolites.ModelSEEDID;  % e.g., Model SEED ID
% Assuming you have two arrays: `reactions` containing reaction formulas and
% `metaboliteList` containing metabolite IDs
% Initialize the metabolite map for indexing
metaboliteMap = containers.Map(metaboliteList, 1:length(metaboliteList));

% Initialize the S matrix
numMetabolites = length(metaboliteList);
numReactions = length(reactions);
S = zeros(numMetabolites, numReactions);

% Initialize the metabolite map for indexing
metaboliteMap = containers.Map(metaboliteList, 1:length(metaboliteList));

% Initialize the S matrix
numMetabolites = length(metaboliteList);
numReactions = length(reactions);
S = zeros(numMetabolites, numReactions);

% Parse each reaction and update S
for i = 1:numReactions
    reaction = reactions{i};
    
    % Identify if the reaction is reversible
    reversible = contains(reaction, '<==>');

    % Split the reaction into substrates and products
    tokens = regexp(reaction, ' <==> | --> ', 'split');
    substrates = strsplit(strtrim(tokens{1}), ' + ');
    products = {};
    if length(tokens) == 2
        products = strsplit(strtrim(tokens{2}), ' + ');
    end

    % Process substrates and products
    for j = 1:length(substrates)
        parts = strsplit(substrates{j}, ' ');
        coeff = -1;  % Default coefficient for substrates
        if length(parts) == 2
            coeff = coeff * str2double(parts{1});
            met_id = parts{2};
        else
            met_id = parts{1};
        end
        if isKey(metaboliteMap, met_id)
            S(metaboliteMap(met_id), i) = S(metaboliteMap(met_id), i) + coeff;
        end
    end
    for j = 1:length(products)
        parts = strsplit(products{j}, ' ');
        coeff = 1;  % Default coefficient for products
        if length(parts) == 2
            coeff = coeff * str2double(parts{1});
            met_id = parts{2};
        else
            met_id = parts{1};
        end
        if isKey(metaboliteMap, met_id)
            S(metaboliteMap(met_id), i) = S(metaboliteMap(met_id), i) + coeff;
        end
    end
end

% Optionally convert S to sparse format
S = sparse(S);

% visualize S matrix
spy(S)

%% check reversibility 
% Assuming 'reactions' is a cell array with your reaction formulas
% You should replace this with the actual import command to load your data
reactions = readtable('GSM.xlsx',sheet = 'iCAC802 Reactions');
reactions = reactions.SEED_FORMULA;
% e.g., reactions = {'cpd00035 + cpd00040 <==> cpd00020 + cpd00033', 'cpd00020 + cpd00054 --> cpd00035 + cpd00145'}; 

% Initialize a vector to store the reversibility information
reversibility = zeros(length(reactions), 1);

% Loop through each reaction formula
for i = 1:length(reactions)
    if contains(reactions{i}, '<==>')
        reversibility(i) = 1; % The reaction is reversible
    elseif contains(reactions{i}, '-->')
        reversibility(i) = 0; % The reaction is irreversible
    end
end

% Now 'reversibility' contains 1s for reversible reactions and 0s for irreversible ones.

%% Create a structured array for your model

reactions = readtable('GSM.xlsx',sheet = 'iCAC802 Reactions');
metabolites = readtable('GSM.xlsx',sheet = 'iCAC802 Metabolites');
% include additional information or rename the item name
model = struct();
model.rxns = reactions.ReactionID;
model.rxnNames = reactions.NAME;
model.subSystems = reactions.SUBSYSTEM;
model.grRules = reactions.GPR;
model.lb = reactions.LowerBound;
model.ub = reactions.UpperBound;
model.ModelSEEDID = metabolites.ModelSEEDID;
model.BiGGID = metabolites.BiGGID;
model.S = S;
model.mets = metabolites.ModelSEEDID;
model.rev = reversibility;

% replace the file name
save('cacmodel.mat', 'model');