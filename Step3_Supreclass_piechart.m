clc
clear
% close all

%% Read significant metabolites

T = readtable('E:\24h brain state\Metabolomics\filtered metabolites0712\MS2 score 0.3\filtered_metabolites_imputed0915.xlsx');

%% Assign missing superclass to "Others"

T.SuperClass = string(T.SuperClass);

idx_missing = ismissing(T.SuperClass) | strlength(strtrim(T.SuperClass))==0;

T.SuperClass(idx_missing) = "Others";

%% Count each superclass

[classNames,~,ic] = unique(T.SuperClass);

counts = accumarray(ic,1);

totalN = sum(counts);

percentages = counts/totalN*100;

%% Merge classes <4% into Others

keep_idx = percentages >= 5 | classNames=="Others";

major_names = classNames(keep_idx);
major_counts = counts(keep_idx);

others_count = sum(counts(~keep_idx));

if others_count > 0

    idxOthers = find(major_names=="Others");

    if isempty(idxOthers)

        major_names(end+1) = "Others";
        major_counts(end+1) = others_count;

    else

        major_counts(idxOthers) = ...
            major_counts(idxOthers) + others_count;

    end

end

%% Recalculate percentages

major_pct = major_counts/sum(major_counts)*100;

%% Sort descending

[major_pct,idx] = sort(major_pct,'descend');

major_names = major_names(idx);
major_counts = major_counts(idx);

%% Output table

ResultTable = table(...
    major_names,...
    major_counts,...
    major_pct,...
    'VariableNames',...
    {'SuperClass','Count','Percentage'});

disp(ResultTable)

%% Save

% Save as Excel file
outputFile = 'SuperClass_Composition.xlsx';
writetable(ResultTable, outputFile);

fprintf('ResultTable has been saved to:\n%s\n', ...
    fullfile(pwd, outputFile));


%% Add percentage labels

txt = compose('%.1f%%',major_pct);

figure

p = pie(major_pct);

text_idx = find(arrayfun(@(x) isa(x,'matlab.graphics.primitive.Text'),p));

for k = 1:length(text_idx)
    p(text_idx(k)).String = txt{k};
end

legend(major_names,...
    'Location','eastoutside')

title('Superclass Composition')



desired_order = [ ...
    "Lipids and lipid-like molecules"
    "Organic acids and derivatives"
    "Organoheterocyclic compounds"
    "Others"
    "Organic oxygen compounds"
    "Nucleosides, nucleotides, and analogues"
    "Benzenoids"
    ];

%% Reorder samples according to the specified order

[~, order_idx] = ismember(desired_order, string(major_names));


if any(order_idx == 0)
    error('Some superclass names were not found in major_names.');
end

major_names_ordered = major_names(order_idx);
major_pct_ordered   = major_pct(order_idx);


txt = compose('%.1f%%', major_pct_ordered);


figure('Color','w');

p = pie(major_pct_ordered);


text_idx = find(arrayfun(@(x) ...
    isa(x,'matlab.graphics.primitive.Text'), p));

for k = 1:length(text_idx)
    p(text_idx(k)).String = txt{k};
end

%% Legend

legend(major_names_ordered, ...
    'Location','eastoutside');

title('Superclass Composition');