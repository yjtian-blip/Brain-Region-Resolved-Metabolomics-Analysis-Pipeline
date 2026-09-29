%%  Heatmap

% Panel 1:
%   Cortex-specific significant metabolites

% Panel 2:
%   Midbrain + Pons-specific significant metabolites
% Panel 3:
%  Commonly altered metabolites across 4 brain regions

% Each row = one metabolite
% Each column = one biological sample

% Data = imputed metabolomics data
% Transformation = log2(Peak Area + 1)
% Visualization = row-wise Z-score

% Sorting:
%   Panel 1:
%       Sort by CTX Log2FC, descending

%   Panel 2:
%       For each metabolite, compare:
%           Midbrain Log2FC
%           Pons Log2FC
%
%       Use the larger of the two Log2FC values
%       as the sorting value
%
%       Then sort descending

clc;
clear;
close all;

%%  Define regions
regions = {'CTX','Hyp','Mid','Pons'};


%%  Read significant-metabolite lists


RegionMet = cell(1,4);

for i = 1:length(regions)

    filename = sprintf('%s_Significant_metabolites.xlsx', ...
        regions{i});

    T = readtable(filename, ...
        'VariableNamingRule','preserve');

    mets = string(T.Metabolites);

    mets = mets(~ismissing(mets));

    RegionMet{i} = unique(mets);

end


% RegionMet:
% RegionMet{1} = CTX significant metabolites
% RegionMet{2} = Hyp significant metabolites
% RegionMet{3} = Mid significant metabolites
% RegionMet{4} = Pons significant metabolites


%%  Create metabolite × region significance matrix


AllMet = unique(vertcat(RegionMet{:}));

MetRegion = false(length(AllMet),4);

for i = 1:4

    MetRegion(:,i) = ...
        ismember(AllMet,RegionMet{i});

end



%%  Identify Cortex-specific metabolites


CTX_specific_idx = ...
    MetRegion(:,1) & ...
    ~MetRegion(:,2) & ...
    ~MetRegion(:,3) & ...
    ~MetRegion(:,4);

CTX_specific_names = ...
    AllMet(CTX_specific_idx);

fprintf('\n============================================\n');

fprintf('Cortex-specific significant metabolites: %d\n', ...
    length(CTX_specific_names));


%%  Identify Midbrain + Pons-specific metabolites



MidPons_specific_idx = ...
    ~MetRegion(:,1) & ...
    ~MetRegion(:,2) & ...
    MetRegion(:,3) & ...
    MetRegion(:,4);

MidPons_specific_names = ...
    AllMet(MidPons_specific_idx);

fprintf('Midbrain + Pons-specific significant metabolites: %d\n', ...
    length(MidPons_specific_names));

fprintf('============================================\n');


%%  Read imputed metabolomics data


raw_file = ...
    'E:\24h brain state\Metabolomics\filtered metabolites0712\MS2 score 0.7\filtered_metabolites_0914\filtered_metabolites_imputed0914.xlsx';

Raw = readtable(raw_file, ...
    'VariableNamingRule','preserve');

Raw_Metabolites = string(Raw.MS2_name);

% Sample columns start from column 9
SampleNames = ...
    string(Raw.Properties.VariableNames(9:end));


%%  Identify samples from each brain region

% ---------------- CTX ----------------

CTX_Ctrl = find( ...
    contains(SampleNames,"CTX") & ...
    contains(lower(SampleNames),"ctr"));

CTX_SD = find( ...
    contains(SampleNames,"CTX") & ...
    contains(lower(SampleNames),"sd"));


% ---------------- Hyp ----------------

Hyp_Ctrl = find( ...
    contains(SampleNames,"Hyp") & ...
    contains(lower(SampleNames),"ctr"));

Hyp_SD = find( ...
    contains(SampleNames,"Hyp") & ...
    contains(lower(SampleNames),"sd"));


% ---------------- Mid ----------------

Mid_Ctrl = find( ...
    contains(SampleNames,"Mid") & ...
    contains(lower(SampleNames),"ctr"));

Mid_SD = find( ...
    contains(SampleNames,"Mid") & ...
    contains(lower(SampleNames),"sd"));


% ---------------- Pons ----------------

Pons_Ctrl = find( ...
    (contains(SampleNames,"Pons") | ...
     contains(SampleNames,"PA")) & ...
    contains(lower(SampleNames),"ctr"));

Pons_SD = find( ...
    (contains(SampleNames,"Pons") | ...
     contains(SampleNames,"PA")) & ...
    contains(lower(SampleNames),"sd"));


%%  Check sample numbers


fprintf('\n========== SAMPLE CHECK ==========\n');

fprintf('CTX Ctrl : %d\n',length(CTX_Ctrl));
fprintf('CTX SD   : %d\n',length(CTX_SD));

fprintf('Hyp Ctrl : %d\n',length(Hyp_Ctrl));
fprintf('Hyp SD   : %d\n',length(Hyp_SD));

fprintf('Mid Ctrl : %d\n',length(Mid_Ctrl));
fprintf('Mid SD   : %d\n',length(Mid_SD));

fprintf('Pons Ctrl: %d\n',length(Pons_Ctrl));
fprintf('Pons SD  : %d\n',length(Pons_SD));

fprintf('==================================\n');


%%  Define FINAL sample order

CTX_Order = [CTX_Ctrl CTX_SD];

Hyp_Order = [Hyp_Ctrl Hyp_SD];

Mid_Order = [Mid_Ctrl Mid_SD];

Pons_Order = [Pons_Ctrl Pons_SD];


All_Order = [ ...
    CTX_Order ...
    Hyp_Order ...
    Mid_Order ...
    Pons_Order];


SampleNames_All = ...
    SampleNames(All_Order);


%% Extract all sample data


% Match all requested metabolites to raw data

[tf_all,loc_all] = ...
    ismember(AllMet,Raw_Metabolites);

AllMet_final = ...
    AllMet(tf_all);

loc_all = ...
    loc_all(tf_all);


% All brain regions × all samples

HeatData_All = ...
    Raw{loc_all,9:end};

HeatData_All = ...
    HeatData_All(:,All_Order);


%%  Extract Cortex-specific data

[tf_CTX,loc_CTX] = ...
    ismember(CTX_specific_names,Raw_Metabolites);

CTX_specific_names_final = ...
    CTX_specific_names(tf_CTX);

loc_CTX = ...
    loc_CTX(tf_CTX);


HeatData_CTX_specific = ...
    Raw{loc_CTX,9:end};

HeatData_CTX_specific = ...
    HeatData_CTX_specific(:,All_Order);


fprintf('\nCortex-specific metabolites found in raw data: %d\n', ...
    size(HeatData_CTX_specific,1));


%%  Extract Midbrain + Pons-specific data


[tf_MidPons,loc_MidPons] = ...
    ismember(MidPons_specific_names,Raw_Metabolites);

MidPons_specific_names_final = ...
    MidPons_specific_names(tf_MidPons);

loc_MidPons = ...
    loc_MidPons(tf_MidPons);


HeatData_MidPons_specific = ...
    Raw{loc_MidPons,9:end};

HeatData_MidPons_specific = ...
    HeatData_MidPons_specific(:,All_Order);


fprintf('Mid + Pons-specific metabolites found in raw data: %d\n', ...
    size(HeatData_MidPons_specific,1));


%%  Calculate Log2FC for sorting
% ---------- Cortex Log2FC ----------

nCTXctrl = length(CTX_Ctrl);

nCTXsd = length(CTX_SD);


CTX_ctrl_data = ...
    HeatData_CTX_specific(:, ...
    1:nCTXctrl);

CTX_SD_data = ...
    HeatData_CTX_specific(:, ...
    nCTXctrl+1:nCTXctrl+nCTXsd);


CTX_Log2FC = ...
    log2( ...
    (mean(CTX_SD_data,2,'omitnan') + eps) ./ ...
    (mean(CTX_ctrl_data,2,'omitnan') + eps));


%%  Sort Cortex-specific metabolites


[CTX_Log2FC_sorted,idx_CTX_sort] = ...
    sort(CTX_Log2FC,'descend');


CTX_specific_names_sorted = ...
    CTX_specific_names_final(idx_CTX_sort);


HeatData_CTX_specific_sorted = ...
    HeatData_CTX_specific(idx_CTX_sort,:);


%%  Calculate Midbrain Log2FC

CTX_total = length(CTX_Order);

Hyp_total = length(Hyp_Order);


Mid_start = ...
    CTX_total + Hyp_total + 1;


Mid_ctrl_start = ...
    Mid_start;

Mid_ctrl_end = ...
    Mid_ctrl_start + length(Mid_Ctrl) - 1;


Mid_SD_start = ...
    Mid_ctrl_end + 1;

Mid_SD_end = ...
    Mid_SD_start + length(Mid_SD) - 1;


Mid_ctrl_data = ...
    HeatData_MidPons_specific(:, ...
    Mid_ctrl_start:Mid_ctrl_end);

Mid_SD_data = ...
    HeatData_MidPons_specific(:, ...
    Mid_SD_start:Mid_SD_end);


Mid_Log2FC = ...
    log2( ...
    (mean(Mid_SD_data,2,'omitnan') + eps) ./ ...
    (mean(Mid_ctrl_data,2,'omitnan') + eps));


%% Calculate Pons Log2FC

Pons_start = ...
    Mid_SD_end + 1;


Pons_ctrl_start = ...
    Pons_start;

Pons_ctrl_end = ...
    Pons_ctrl_start + length(Pons_Ctrl) - 1;


Pons_SD_start = ...
    Pons_ctrl_end + 1;

Pons_SD_end = ...
    Pons_SD_start + length(Pons_SD) - 1;


Pons_ctrl_data = ...
    HeatData_MidPons_specific(:, ...
    Pons_ctrl_start:Pons_ctrl_end);

Pons_SD_data = ...
    HeatData_MidPons_specific(:, ...
    Pons_SD_start:Pons_SD_end);


Pons_Log2FC = ...
    log2( ...
    (mean(Pons_SD_data,2,'omitnan') + eps) ./ ...
    (mean(Pons_ctrl_data,2,'omitnan') + eps));


%%  Determine sorting Log2FC for Midbrain + Pons


MidPons_sort_Log2FC = ...
    max([Mid_Log2FC Pons_Log2FC],[],2);


%%  Sort Midbrain + Pons-specific metabolites


[MidPons_sort_Log2FC_sorted, ...
    idx_MidPons_sort] = ...
    sort(MidPons_sort_Log2FC,'descend');


%% Reorder metabolite names

MidPons_specific_names_sorted = ...
    MidPons_specific_names_final(idx_MidPons_sort);


%% Reorder heatmap data

HeatData_MidPons_specific_sorted = ...
    HeatData_MidPons_specific(idx_MidPons_sort,:);


%% Reorder individual Log2FC values

Mid_Log2FC_sorted = ...
    Mid_Log2FC(idx_MidPons_sort);

Pons_Log2FC_sorted = ...
    Pons_Log2FC(idx_MidPons_sort);


%%  Log2 transformation

CTX_log = ...
    log2(HeatData_CTX_specific_sorted + 1);


MidPons_log = ...
    log2(HeatData_MidPons_specific_sorted + 1);


%%  Row-wise Z-score

CTX_z = ...
    zscore(CTX_log,0,2);

MidPons_z = ...
    zscore(MidPons_log,0,2);


%%  Remove rows with NaN / Inf Z-score

% Keep metabolite names and Log2FC values aligned

valid_CTX = ...
    all(isfinite(CTX_z),2);

valid_MidPons = ...
    all(isfinite(MidPons_z),2);


%% ---------- Cortex ----------

CTX_z = ...
    CTX_z(valid_CTX,:);

CTX_specific_names_sorted = ...
    CTX_specific_names_sorted(valid_CTX);

CTX_Log2FC_sorted = ...
    CTX_Log2FC_sorted(valid_CTX);


%% ---------- Mid + Pons ----------

MidPons_z = ...
    MidPons_z(valid_MidPons,:);

MidPons_specific_names_sorted = ...
    MidPons_specific_names_sorted(valid_MidPons);

MidPons_sort_Log2FC_sorted = ...
    MidPons_sort_Log2FC_sorted(valid_MidPons);

Mid_Log2FC_sorted = ...
    Mid_Log2FC_sorted(valid_MidPons);

Pons_Log2FC_sorted = ...
    Pons_Log2FC_sorted(valid_MidPons);


%%  Create brain-region boundary positions

nCTX  = length(CTX_Order);

nHyp  = length(Hyp_Order);

nMid  = length(Mid_Order);

nPons = length(Pons_Order);


%%  Plot two-panel heatmap


figure( ...
    'Position',[50 50 1800 1000], ...
    'Color','w');


tiledlayout(2,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% PANEL 1 Cortex-specific metabolites

ax1 = nexttile;


imagesc(ax1,CTX_z);


colormap(ax1,redbluecmap);

caxis(ax1,[-2 2]);


%% ---------- Y axis ----------

yticks(ax1, ...
    1:length(CTX_specific_names_sorted));

yticklabels(ax1, ...
    CTX_specific_names_sorted);


%% ---------- X axis ----------

xticks(ax1, ...
    1:length(SampleNames_All));

xticklabels(ax1, ...
    SampleNames_All);

xtickangle(ax1,90);


xlabel(ax1,'Samples');

ylabel(ax1,'Metabolites');


title(ax1, ...
    sprintf( ...
    'Cortex-specific significant metabolites (n=%d)', ...
    length(CTX_specific_names_sorted)), ...
    'FontSize',14, ...
    'FontWeight','bold');


set(ax1, ...
    'FontSize',9, ...
    'LineWidth',1, ...
    'Box','off');


%% ---------- Brain-region separators ----------

hold(ax1,'on');


xline(ax1, ...
    nCTX+0.5, ...
    'k-', ...
    'LineWidth',1.5);


xline(ax1, ...
    nCTX+nHyp+0.5, ...
    'k-', ...
    'LineWidth',1.5);


xline(ax1, ...
    nCTX+nHyp+nMid+0.5, ...
    'k-', ...
    'LineWidth',1.5);


hold(ax1,'off');


%% PANEL 2 Midbrain + Pons-specific metabolites

ax2 = nexttile;


imagesc(ax2,MidPons_z);


colormap(ax2,redbluecmap);

caxis(ax2,[-2 2]);


%% ---------- Y axis ----------

yticks(ax2, ...
    1:length(MidPons_specific_names_sorted));

yticklabels(ax2, ...
    MidPons_specific_names_sorted);


%% ---------- X axis ----------

xticks(ax2, ...
    1:length(SampleNames_All));

xticklabels(ax2, ...
    SampleNames_All);

xtickangle(ax2,90);


xlabel(ax2,'Samples');

ylabel(ax2,'Metabolites');


title(ax2, ...
    sprintf( ...
    'Midbrain + Pons-specific significant metabolites (n=%d)', ...
    length(MidPons_specific_names_sorted)), ...
    'FontSize',14, ...
    'FontWeight','bold');


set(ax2, ...
    'FontSize',9, ...
    'LineWidth',1, ...
    'Box','off');


%% ---------- Brain-region separators ----------

hold(ax2,'on');


xline(ax2, ...
    nCTX+0.5, ...
    'k-', ...
    'LineWidth',1.5);


xline(ax2, ...
    nCTX+nHyp+0.5, ...
    'k-', ...
    'LineWidth',1.5);


xline(ax2, ...
    nCTX+nHyp+nMid+0.5, ...
    'k-', ...
    'LineWidth',1.5);


hold(ax2,'off');




cb = colorbar(ax2);

cb.Label.String = ...
    'Row Z-score';

cb.FontSize = 11;



fprintf('\n============================================\n');

fprintf('FINAL SUMMARY\n');

fprintf('--------------------------------------------\n');

fprintf('Cortex-specific metabolites : %d\n', ...
    length(CTX_specific_names_sorted));

fprintf('Mid + Pons-specific metabolites : %d\n', ...
    length(MidPons_specific_names_sorted));

fprintf('--------------------------------------------\n');

fprintf('Panel 1 sorting:\n');

fprintf('  CTX Log2FC descending\n');

fprintf('\nPanel 2 sorting:\n');

fprintf('  max(Mid Log2FC, Pons Log2FC) descending\n');

fprintf('--------------------------------------------\n');

fprintf('Cortex sample columns:\n');

fprintf('  CTX  Ctrl = %d\n',length(CTX_Ctrl));

fprintf('  CTX  SD   = %d\n',length(CTX_SD));


fprintf('Hyp sample columns:\n');

fprintf('  Hyp  Ctrl = %d\n',length(Hyp_Ctrl));

fprintf('  Hyp  SD   = %d\n',length(Hyp_SD));


fprintf('Mid sample columns:\n');

fprintf('  Mid  Ctrl = %d\n',length(Mid_Ctrl));

fprintf('  Mid  SD   = %d\n',length(Mid_SD));


fprintf('Pons sample columns:\n');

fprintf('  Pons Ctrl = %d\n',length(Pons_Ctrl));

fprintf('  Pons SD   = %d\n',length(Pons_SD));


fprintf('============================================\n');

%% Heatmap for commonly-altered metabolites 
regions = {'Hyp','Mid','Pons','CTX'};

%  Step 1. Read each region and store metabolite + Log2FC

RegionTables = struct();

for i = 1:length(regions)

    region = regions{i};

    filename = sprintf('%s_Significant_Metabolites.xlsx', region);

    T = readtable(filename, 'VariableNamingRule','preserve');

    % Store the original table
    RegionTables.(region) = T;

end


% Step 2. Collect all unique metabolites


AllMet = strings(0,1);

for i = 1:length(regions)

    region = regions{i};
    T = RegionTables.(region);

    AllMet = [AllMet; string(T.Metabolites)];

end

UniqueMet = unique(AllMet);


%  Step 3. Create output table

Result = table(UniqueMet, ...
    'VariableNames', {'Metabolites'});


% Step 4. Add Log2FC from each brain region

for i = 1:length(regions)

    region = regions{i};
    T = RegionTables.(region);

    MetNames = string(T.Metabolites);

    % Find matching metabolites
    [tf, loc] = ismember(UniqueMet, MetNames);

    % Initialize with NaN
    Log2FC = nan(height(Result),1);

    % Fill matched metabolites
    Log2FC(tf) = T.Log2FC(loc(tf));

    % Add to output table
    Result.(sprintf('%s_Log2FC',region)) = Log2FC;

end


% Step 5. Determine whether metabolite is present in all 4 regions


Shared4 = ...
    ~isnan(Result.Hyp_Log2FC) & ...
    ~isnan(Result.Mid_Log2FC) & ...
    ~isnan(Result.Pons_Log2FC) & ...
    ~isnan(Result.CTX_Log2FC);


% Add logical column LAST
Result.Shared4 = Shared4;


% Step 6. Keep only metabolites shared by all 4 regions

Shared4_Result = Result(Result.Shared4,:);


%%  Save


writetable(Shared4_Result, ...
    'Shared_Significant_Metabolites_4Regions_0922.xlsx');

%%  Extract significantly altered metabolites across 4 brain regions
shared_file = 'E:\24h brain state\Metabolomics\filtered metabolites0712\MS2 score 0.7\sigdata\Shared_Significant_Metabolites_4Regions_0922.xlsx';
SharedT = readtable(shared_file, 'VariableNamingRule', 'preserve');

SharedMetNames = string(SharedT.Metabolites); 

%% Read imputed data
raw_file = 'E:\24h brain state\Metabolomics\filtered metabolites0712\MS2 score 0.7\filtered_metabolites_0914\filtered_metabolites_imputed0914.xlsx';
Raw = readtable(raw_file, 'VariableNamingRule', 'preserve');

Raw_Metabolites = string(Raw. MS2_name);

[tf, loc] = ismember(SharedMetNames, Raw_Metabolites);

matched_loc = loc(tf);
Final_Met_Names = SharedMetNames(tf);


HeatData = Raw{matched_loc, 9:end};

SampleNames = Raw.Properties.VariableNames(9:end);

%% Data normalization

HeatData_log = log2(HeatData + 1);


HeatData_z = zscore(HeatData_log, 0, 2); 

%%  Save results


Result_Zscore = array2table(HeatData_z, 'VariableNames', SampleNames);

%
Result_Zscore = addvars(Result_Zscore, Final_Met_Names, 'Before', 1, 'NewVariableNames', 'MS2Name');

save_zscore_filename = 'Shared_Metabolites_Zscore_Data2.xlsx';
writetable(Result_Zscore, save_zscore_filename);

fprintf('==== 处理完成 ====\n');
fprintf('成功提取并标准化了 %d 个共享代谢物的 Peak Area。\n', length(Final_Met_Names));
fprintf('Z-score 矩阵已保存至: %s\n', save_zscore_filename);

SampleNames = string(SampleNames);

%% Define brain regions
CTX = contains(SampleNames,"CTX");
Hyp = contains(SampleNames,"Hyp");
Mid = contains(SampleNames,"Mid");
PA  = contains(SampleNames,"PA");

%% Define Control and SD samples
CTX_Ctrl = find(contains(SampleNames,"CTX") & ...
                contains(lower(SampleNames),"ctr"));

CTX_SD   = find(contains(SampleNames,"CTX") & ...
                contains(lower(SampleNames),"sd"));

Hyp_Ctrl = find(contains(SampleNames,"Hyp") & ...
                contains(lower(SampleNames),"ctr"));

Hyp_SD   = find(contains(SampleNames,"Hyp") & ...
                contains(lower(SampleNames),"sd"));

Mid_Ctrl = find(contains(SampleNames,"Mid") & ...
                contains(lower(SampleNames),"ctr"));

Mid_SD   = find(contains(SampleNames,"Mid") & ...
                contains(lower(SampleNames),"sd"));

PA_Ctrl  = find(contains(SampleNames,"PA") & ...
                contains(lower(SampleNames),"ctr"));

PA_SD    = find(contains(SampleNames,"PA") & ...
                contains(lower(SampleNames),"sd"));

%% Reorder samples
newOrder = [...
    CTX_Ctrl CTX_SD ...
    Hyp_Ctrl Hyp_SD ...
    Mid_Ctrl Mid_SD ...
    PA_Ctrl  PA_SD];

%% Keep original data for Log2FC calculation
HeatData_original = HeatData;

%% Calculate Log2FC for each brain region
CTX_Log2FC = log2( ...
    (mean(HeatData_original(:,CTX_SD),2,'omitnan') + eps) ./ ...
    (mean(HeatData_original(:,CTX_Ctrl),2,'omitnan') + eps));

Hyp_Log2FC = log2( ...
    (mean(HeatData_original(:,Hyp_SD),2,'omitnan') + eps) ./ ...
    (mean(HeatData_original(:,Hyp_Ctrl),2,'omitnan') + eps));

Mid_Log2FC = log2( ...
    (mean(HeatData_original(:,Mid_SD),2,'omitnan') + eps) ./ ...
    (mean(HeatData_original(:,Mid_Ctrl),2,'omitnan') + eps));

PA_Log2FC = log2( ...
    (mean(HeatData_original(:,PA_SD),2,'omitnan') + eps) ./ ...
    (mean(HeatData_original(:,PA_Ctrl),2,'omitnan') + eps));

%% Calculate maximum Log2FC across the four brain regions
Max_Log2FC = max([...
    CTX_Log2FC ...
    Hyp_Log2FC ...
    Mid_Log2FC ...
    PA_Log2FC], [], 2, 'omitnan');

%% Sort metabolites according to Max(Log2FC)
[Max_Log2FC_sorted, sort_idx] = sort(Max_Log2FC, 'descend');

%% Reorder heatmap data
HeatData_z_ordered = HeatData_z(:,newOrder);

HeatData_z_final = HeatData_z_ordered(sort_idx,:);

Final_Met_Names_sorted = Final_Met_Names(sort_idx);

%% Reorder sample names
SampleNames_ordered = SampleNames(newOrder);

%% Plot heatmap
figure('Position',[150 150 1200 800]);

imagesc(HeatData_z_final);

colormap(redbluecmap);

caxis([-2 2]);

colorbar;

%% Y axis
yticks(1:length(Final_Met_Names_sorted));
yticklabels(Final_Met_Names_sorted);

%% X axis
xticks(1:length(SampleNames_ordered));
xticklabels(SampleNames_ordered);
xtickangle(90);

ylabel('Metabolites');
xlabel('Samples');

title('Z-score normalized metabolite abundance');
%% Create sorting check table

SortCheck = table( ...
    Final_Met_Names(:), ...
    CTX_Log2FC(:), ...
    Hyp_Log2FC(:), ...
    Mid_Log2FC(:), ...
    PA_Log2FC(:), ...
    Max_Log2FC(:), ...
    'VariableNames', { ...
    'Metabolite', ...
    'CTX_Log2FC', ...
    'Hyp_Log2FC', ...
    'Mid_Log2FC', ...
    'PA_Log2FC', ...
    'Max_Log2FC'});

SortCheck = SortCheck(sort_idx,:);

disp(SortCheck);