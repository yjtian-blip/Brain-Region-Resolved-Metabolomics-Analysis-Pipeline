clc
clear
close all

%% Read significant metabolites

T = readtable('G:\20260715\filtered_metabolites_imputed.xlsx');

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

keep_idx = percentages >= 4 | classNames=="Others";

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

writetable(ResultTable,...
    'Superclass_Composition.xlsx');


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