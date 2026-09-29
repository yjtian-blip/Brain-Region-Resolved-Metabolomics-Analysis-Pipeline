clc;
clear;
% close all;

%% Read table

T = readtable('E:\24h brain state\Metabolomics\filtered metabolites0712\MS2 score 0.7\CTX\pathway_results.csv');
Pathway = string(T.Pathway);
Impact  = T.Impact;
AdjustP = T.AdjustP;
%Sort metabolites by AdjustP
[AdjustP, idx] = sort(AdjustP, 'descend');
Pathway = Pathway(idx);
Impact  = Impact(idx);

topN = min(7, length(Pathway)); 
Pathway_top = Pathway(1:topN);
Impact_top  = Impact(1:topN);
AdjustP_top = AdjustP(1:topN);


Pathway_top = flip(Pathway_top);
Impact_top  = flip(Impact_top);
AdjustP_top = flip(AdjustP_top);

Y = 1:topN;



BubbleSize = 50 + 800 * Impact_top; 



figure('Position',[200 100 800 450]);


scatter(AdjustP_top,...
        Y,...
        BubbleSize,...
        AdjustP_top,...
        'filled',...
        'MarkerEdgeColor','k');
%% Set axis

yticks(Y);
yticklabels(Pathway_top);
xlabel('AdjustP'); 
ylabel('');
set(gca,...
    'FontSize',11,...
    'LineWidth',1.2,...
    'TickLength',[0 0],...
    'YLim',[0.5 topN+0.5]); 

R = linspace(200/255, 200/255, 256)';   
G = linspace(200/255,  12/255, 256)';   
B = linspace(200/255,  12/255, 256)';  
custom_cmap = [R, G, B];
colormap(custom_cmap);

cb = colorbar;
cb.Label.String = 'Adjusted P';

if max(AdjustP_top) > min(AdjustP_top)
    caxis([min(AdjustP_top) max(AdjustP_top)]); 
  
end


grid on;
box off;
set(gca,'Color',[0.95 0.95 0.95]);
title('Top 5 KEGG Pathway Enrichment');
