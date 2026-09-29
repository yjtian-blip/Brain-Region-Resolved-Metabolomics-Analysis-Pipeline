% UpSet plot for differential metabolites across 4 brain regions
%Note: Compiled significant metabolites across all four brain regions into total_sig_metabolites.xlsx

clear;
clc;
close all;


%% Read differential metabolite lists


file_path = 'E:\24h brain state\Metabolomics\filtered metabolites0712\MS2 score 0.7\sigdata\total_sig_metabolites.xlsx';

T = readtable(file_path);


CTX  = string(T{:,1});
Hyp  = string(T{:,2});
Mid  = string(T{:,3});
Pons = string(T{:,4});


%% Remove missing / empty cells

CTX  = CTX(~ismissing(CTX)  & strlength(strtrim(CTX)) > 0);
Hyp  = Hyp(~ismissing(Hyp)  & strlength(strtrim(Hyp)) > 0);
Mid  = Mid(~ismissing(Mid)  & strlength(strtrim(Mid)) > 0);
Pons = Pons(~ismissing(Pons) & strlength(strtrim(Pons)) > 0);


%%  Remove duplicated metabolites

CTX  = unique(CTX);
Hyp  = unique(Hyp);
Mid  = unique(Mid);
Pons = unique(Pons);


%%  Set size

region_names = ["Cortex","Hypothalamus","Midbrain","Pons"];

set_size = [
    length(CTX), ...
    length(Hyp), ...
    length(Mid), ...
    length(Pons)
];


%% Display

fprintf('\n====================================\n');
fprintf('Differential metabolite set size\n');
fprintf('====================================\n');

for i = 1:4
    fprintf('%s: %d\n', ...
        region_names(i),set_size(i));
end


%% Define all unique metabolites

all_metabolites = unique([
    CTX;
    Hyp;
    Mid;
    Pons
]);


n_metabolites = length(all_metabolites);


%% Construct membership matrix
%
% Columns:
%       CTX | Hyp | Mid | Pons
%
% Each row = one metabolite
%
% 1 = present
% 0 = absent
% ==========================================================

membership = false(n_metabolites,4);

membership(:,1) = ismember(all_metabolites,CTX);
membership(:,2) = ismember(all_metabolites,Hyp);
membership(:,3) = ismember(all_metabolites,Mid);
membership(:,4) = ismember(all_metabolites,Pons);


%%  Calculate SET SIZE

% Total number of differential metabolites in each region


set_size = sum(membership,1);

region_names = ["CTX","Hyp","Mid","Pons"];


fprintf('\n============================================\n');
fprintf('Set size\n');
fprintf('============================================\n');

for i = 1:4
    fprintf('%s: %d differential metabolites\n', ...
        region_names(i),set_size(i));
end

combination_matrix = dec2bin(1:15) - '0';

combination_matrix = logical(combination_matrix);


%%  Calculate EXACT intersection size


n_combinations = size(combination_matrix,1);

intersection_size = zeros(n_combinations,1);

for i = 1:n_combinations

    combination = combination_matrix(i,:);

    % Metabolites that have exactly this membership pattern
    exact_idx = all(membership == combination,2);

    intersection_size(i) = sum(exact_idx);

end


%%  Remove empty intersections


nonzero_idx = intersection_size > 0;

intersection_size = intersection_size(nonzero_idx);
combination_matrix = combination_matrix(nonzero_idx,:);


%%  Sort intersections by size


[intersection_size,sort_idx] = ...
    sort(intersection_size,'descend');

combination_matrix = combination_matrix(sort_idx,:);


%% Display intersection results


fprintf('\n============================================\n');
fprintf('Intersection size\n');
fprintf('============================================\n');

for i = 1:length(intersection_size)

    current_combination = combination_matrix(i,:);

    current_regions = region_names(current_combination);

    combination_name = strjoin(current_regions,' + ');

    fprintf('%-25s %d\n', ...
        combination_name, ...
        intersection_size(i));

end


fprintf('\n============================================\n');
fprintf('Check\n');
fprintf('============================================\n');

fprintf('Total unique metabolites = %d\n',n_metabolites);
fprintf('Sum of intersections     = %d\n',sum(intersection_size));


%% Create UpSet plot


figure( ...
    'Color','w', ...
    'Position',[100 100 1200 800]);




% Main intersection bar plot
ax1 = axes( ...
    'Position',[0.32 0.50 0.62 0.42]);

hold(ax1,'on');


% Upper bar plot: Intersection size

x = 1:length(intersection_size);

bar(ax1,x,intersection_size, ...
    'FaceAlpha',0.85, ...
    'EdgeColor','none');

ylabel(ax1,'Intersection size', ...
    'FontSize',14, ...
    'FontWeight','bold');

set(ax1, ...
    'FontSize',12, ...
    'LineWidth',1.2, ...
    'Box','off');

xlim(ax1,[0.4 length(x)+0.6]);

xticks(ax1,[]);

title(ax1, ...
    'Differential metabolite intersections', ...
    'FontSize',16, ...
    'FontWeight','bold');




ax2 = axes( ...
    'Position',[0.32 0.12 0.62 0.27]);

hold(ax2,'on');


% Plot background rows

for r = 1:4

    yline(ax2,r, ...
        'Color',[0.90 0.90 0.90], ...
        'LineWidth',20);

end


%% Plot membership matrix


for i = 1:length(intersection_size)

    active = find(combination_matrix(i,:));
    inactive = find(~combination_matrix(i,:));

    % Draw connecting line when more than one set
    if length(active) > 1

        plot(ax2, ...
            [i i], ...
            [min(active) max(active)], ...
            'k-', ...
            'LineWidth',2);

    end

    % Plot inactive grey dots
    plot(ax2, ...
        repmat(i,1,length(inactive)), ...
        inactive, ...
        'o', ...
        'MarkerSize',9, ...
        'MarkerFaceColor',[0.80 0.80 0.80], ...
        'MarkerEdgeColor','none');

    % Plot active black dots
    plot(ax2, ...
        repmat(i,1,length(active)), ...
        active, ...
        'o', ...
        'MarkerSize',10, ...
        'MarkerFaceColor','k', ...
        'MarkerEdgeColor','k');

end


% Format membership matrix


xlim(ax2,[0.4 length(x)+0.6]);

ylim(ax2,[0.5 4.5]);

yticks(ax2,1:4);

yticklabels(ax2,region_names);

set(ax2, ...
    'YDir','reverse', ...
    'FontSize',13, ...
    'LineWidth',1.2, ...
    'Box','off');

xticks(ax2,[]);

ylabel(ax2,'');

%% LEFT panel: Set size

ax3 = axes('Position',[0.08 0.12 0.20 0.27]);
hold(ax3,'on');

% ---------------------------------------------------------
% Horizontal bars extending from 0 toward the left
% ---------------------------------------------------------
barh(ax3, 1:4, -set_size, ...
    'FaceAlpha',0.85, ...
    'EdgeColor','none');

% ---------------------------------------------------------
% Y axis
% ---------------------------------------------------------
set(ax3, ...
    'YDir','reverse', ...
    'YAxisLocation','right', ...
    'FontSize',13, ...
    'LineWidth',1.2, ...
    'Box','off');

yticks(ax3,1:4);
yticklabels(ax3,region_names);

% ---------------------------------------------------------
% X axis
% ---------------------------------------------------------
max_set = max(set_size);
xmax = ceil(max_set/10)*10;

xlim(ax3,[-140 0]);

xt = -140:20:0;
xticks(ax3,xt);
xticklabels(ax3,string(abs(xt)));

% Rotate tick labels by 45 degrees
xtickangle(ax3,45);

xlabel(ax3,'Set size', ...
    'FontSize',13, ...
    'FontWeight','bold');

% ---------------------------------------------------------
% Add numbers at the left end of bars
% ---------------------------------------------------------
for i = 1:4

    text(ax3, ...
        -set_size(i)-2, ...
        i, ...
        sprintf('%d',set_size(i)), ...
        'HorizontalAlignment','right', ...
        'VerticalAlignment','middle', ...
        'FontSize',11);

end
% Add numbers above intersection bars


max_y = max(intersection_size);

for i = 1:length(intersection_size)

    text(ax1, ...
        i, ...
        intersection_size(i), ...
        sprintf('%d',intersection_size(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontSize',10);

end


%% Final formatting


linkaxes([ax1 ax2],'x');

hold(ax1,'off');
hold(ax2,'off');
hold(ax3,'off');

saveas(gcf,'UpSet_differential_metabolites.png');
