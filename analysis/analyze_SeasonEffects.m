
clear

% add paths
addpath('../lib/cbrewer/')
addpath('../lib/fdr_bh/')

% load data
load('../KEGG/parsed_KO_paths_ABnC_cutoff_02.mat')
load('../data/meta_full.mat', 'meta')


[~,i,j]=intersect(meta.hetero.metagenome_genes,out.koList);

meta.hetero.metagenomes= meta.hetero.metagenomes(i,:);
meta.hetero.metagenome_genes = meta.hetero.metagenome_genes(i);

z = sum(out.K(j,:)) > 0;
P = out.K(j,z);
paths = out.path(z);
X = cell2mat(cellfun(@(x) x./sum(x), num2cell(meta.hetero.metagenomes,1),'uni',0));
meta.hetero.meta_paths_rel = P'*X;
meta.hetero.paths = paths;

idx22 = meta.samples_info.Filter == 0.22;
idx11 = meta.samples_info.Filter == 11;


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
 
 %% see if AIC lowers with inclusion of season information 
 
%   for i = 1:length(meta.hetero.paths)
%      %x = meta.hetero.meta_paths_rel(i,idx11);
%      %sm = find(strcmp(meta.samples_info.Season(idx11),'Sm'));
%      %%wt = find(strcmp(meta.samples_info.Season(idx11),'Wt'));
%      %sp = find(strcmp(meta.samples_info.Season(idx11),'Sp'));
%      y = meta.hetero.meta_paths_rel(i,:)';
%      %p_anova_11(i) = anova1(y,meta.samples_info.Season(:),'off');
%   
%      
%      %Xf = dm.Filter(idx11,:);
%      Xd = meta.samples_info.Depth_m(:);
%      Xs = dm.Season(:,:);
%      Xf = dm.Filter;
%      %X = dm.Season(idx11,:);
%      
%      mdl = fitglm([Xf Xd],y);
%      aic = mdl.ModelCriterion.AICc;
%      
%      mdl2 = fitglm([Xf,Xd,Xs],y);
%      aic2 = mdl2.ModelCriterion.AICc;
%      
%      if aic2 < aic
%          season(i) = true;
%      else
%          season(i) = false;
%      end
%      
%      
%  end
%  season_full = season;
%  
 
 
%  for i = 1:length(meta.hetero.paths)
%      p_anova_22(i) = anova1(meta.hetero.meta_paths_rel(i,idx22),meta.samples_info.Season(idx22),'off');
%      p_anova_11(i) = anova1(meta.hetero.meta_paths_rel(i,idx11),meta.samples_info.Season(idx11),'off');
%  end
 
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
 pa.h = h;
 pa.season = season;
 pa.padj = padj;
 pa.pvals = p_anova_11;
 
 
 
 
 
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
 
 
 fl.h = h;
 fl.season = season;
 fl.padj = padj;
 fl.pvals = p_anova_22;
 
 
 %% PLOT DATA
 z1 = find(pa.h);
 z2 = find(fl.h);
 
 X = meta.hetero.meta_paths_rel(z1,idx11);
 D = squareform(pdist(X));
 Z = linkage(D);
 o = optimalleaforder(Z,D);
 z1 = z1(o);
 
 
 figure();
 b = 1;
 cmap = [cbrewer('qual','Set1',3)];
 padj = pa.padj(z1)';
 for i = 1:length(z1)
    subplot(9,5,b)
    y = meta.hetero.meta_paths_rel(z1(i),idx11)';
    g = meta.samples_info.Season(idx11);
    boxplot(y,g,'plotstyle','compact','colors',cmap,'GroupOrder',{'Wt','Sp','Sm'},'Labels',{'' '' ''},'Widths',1.4);
    b = b+1;
    title(num2str(i));
    set(gca,'XTick',[]);
    text(1.5,min(y)+0.01*min(y), ['{\sl{P}}_{adj} = ', num2str(round(pa.padj(z1(i)),3))])
    %set(gcf,'Color',[0.6 0.6 0.6]);

 end
 
  for i = 1:length(z2)
    subplot(9,5,b)
    y = meta.hetero.meta_paths_rel(z2(i),idx22)';
    g = meta.samples_info.Season(idx22);
    boxplot(y,g,'plotstyle','compact','colors',cmap,'GroupOrder',{'Wt','Sp','Sm'},'Labels',{'' '' ''},'Widths',1.4);
    b = b+1;
    title(num2str(i));
    set(gca,'XTick',[]);
        text(1.5,min(y)+0.01*min(y), ['{\sl{P}}_{adj} = ', num2str(round(fl.padj(z2(i)),3))])

    %set(gcf,'Color',[0.6 0.6 0.6]);

 end
 

 
 