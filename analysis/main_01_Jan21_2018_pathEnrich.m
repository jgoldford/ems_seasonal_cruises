% Script to analysis the abundance of heterotrophic pathways at various depths 
% Author: Joshua E. Goldford


clear

% add paths
addpath('../lib/cbrewer/')
addpath('../lib/fdr_bh/')

% load data
load('../KEGG/parsed_KO_paths_ABnC_cutoff_02.mat')
load('../data/meta_full.mat', 'meta')

% remove genes in the metagenome structure that are not used in KEGG pathways
[~,i,j]=intersect(meta.hetero.metagenome_genes,out.koList);
meta.hetero.metagenomes= meta.hetero.metagenomes(i,:);
meta.hetero.metagenome_genes = meta.hetero.metagenome_genes(i);

% remove pathways with no genes from PiCRUST analysis
z = sum(out.K(j,:)) > 0;
% create a pathway-to-gene correspondence matrix (used in the following step)
P = out.K(j,z);
% list the pathways by ID
paths = out.path(z);

% normalize metagenomes
X = cell2mat(cellfun(@(x) x./sum(x), num2cell(meta.hetero.metagenomes,1),'uni',0));
% for each pathway, sum up all relative abundances of each participating gene, creating a matrix of pathways-to-samples that encode the 
% pathway-level fractional abundances. Append this to the struct containing all data 
meta.hetero.meta_paths_rel = P'*X;
meta.hetero.paths = paths;


% find indicies for filter sizes
idx22 = meta.samples_info.Filter == 0.22;
idx11 = meta.samples_info.Filter == 11;


% for each pathway and filter size, compute the correlation and slopes between pathway-level relative abundances and sample depth
for i = 1:length(meta.hetero.paths)

    % compute pearson correlation between pathway relative abundance and sample depth (r = correlation coefficient, p = p-value)
    [r,p]=corr(meta.hetero.meta_paths_rel(i,idx22)',meta.samples_info.Depth_m(idx22));
    R22(i) = r;
    P22(i) = p;
    M22(i) = mean(meta.hetero.meta_paths_rel(i,idx22));
    % fit a linear model, and compute the slope between pathway relative abundance and sample depth
    B22{i} = glmfit(meta.hetero.meta_paths_rel(i,idx22)',meta.samples_info.Depth_m(idx22));
   
    % compute pearson correlation between pathway relative abundance and sample depth (r = correlation coefficient, p = p-value)
    [r,p]=corr(meta.hetero.meta_paths_rel(i,idx11)',meta.samples_info.Depth_m(idx11));
    R11(i) = r;
    P11(i) = p;
    M11(i) = mean(meta.hetero.meta_paths_rel(i,idx11));
   % fit a linear model, and compute the slope between pathway relative abundance and sample depth
    B11{i} = glmfit(meta.hetero.meta_paths_rel(i,idx11)',meta.samples_info.Depth_m(idx11));

    
end
 


%% analyze seasons 
 [~,~,L]=unique(meta.samples_info.Season);
 dm.Season = zeros(size(meta.hetero.meta_paths_rel,2),max(L));
 for row =1:size(dm.Season,1)
    dm.Season(row,L(row)) = 1; 
 end

 [~,~,L]=unique(meta.samples_info.Filter);
 dm.Filter = zeros(size(meta.hetero.meta_paths_rel,2),max(L));
 for row =1:size(dm.Season,1)
    dm.Filter(row,L(row)) = 1; 
 end
 
 
 for i = 1:length(meta.hetero.paths)
     p_anova_22(i) = anova1(meta.hetero.meta_paths_rel(i,idx22),meta.samples_info.Season(idx22),'off');
     p_anova_11(i) = anova1(meta.hetero.meta_paths_rel(i,idx11),meta.samples_info.Season(idx11),'off');
 end
 
 for i = 1:length(meta.hetero.paths)
     %x = meta.hetero.meta_paths_rel(i,idx11);
     %sm = find(strcmp(meta.samples_info.Season(idx11),'Sm'));
     %%wt = find(strcmp(meta.samples_info.Season(idx11),'Wt'));
     %sp = find(strcmp(meta.samples_info.Season(idx11),'Sp'));
     y = meta.hetero.meta_paths_rel(i,idx11)';
     p_anova_11(i) = anova1(y,meta.samples_info.Season(idx11),'off');
  
     
     %Xf = dm.Filter(idx11,:);
     Xd = meta.samples_info.Depth_m(idx11);
     Xs = dm.Season(idx11,:);
     
     %X = dm.Season(idx11,:);
     
     
     mdl = fitglm(Xd,y);
     aic = mdl.ModelCriterion.AICc;
     
     mdl2 = fitglm([Xd,Xs],y);
     aic2 = mdl2.ModelCriterion.AICc;
     
     if aic2 < aic
         season(i) = true;
     else
         season(i) = false;
     end
     
     
 end
 [h,~,~,padj]=fdr_bh(p_anova_11);
 [~,l]=sort(padj);
 
 
 
 for i = 1:length(meta.hetero.paths)
     %x = meta.hetero.meta_paths_rel(i,idx11);
     %sm = find(strcmp(meta.samples_info.Season(idx11),'Sm'));
     %%wt = find(strcmp(meta.samples_info.Season(idx11),'Wt'));
     %sp = find(strcmp(meta.samples_info.Season(idx11),'Sp'));
     y = meta.hetero.meta_paths_rel(i,idx22)';
     p_anova_22(i) = anova1(y,meta.samples_info.Season(idx22),'off');
  
     
     %Xf = dm.Filter(idx11,:);
     Xd = meta.samples_info.Depth_m(idx22);
     Xs = dm.Season(idx22,:);
     
     %X = dm.Season(idx11,:);
     
     
     mdl = fitglm(Xd,y);
     aic = mdl.ModelCriterion.AICc;
     
     mdl2 = fitglm([Xd,Xs],y);
     aic2 = mdl2.ModelCriterion.AICc;
     
     if aic2 < aic
         season(i) = true;
     else
         season(i) = false;
     end
     
     
 end
 [h,~,~,padj]=fdr_bh(p_anova_22);
 [~,l]=sort(padj);
 
 
 
