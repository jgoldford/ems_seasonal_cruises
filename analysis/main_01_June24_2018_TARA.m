% analysis of TARA oceans data (at the pathway level)
% June-23-2018
% Joshua Goldford

clear
% load the data
load('../TARA/tara_data_parsed.mat')

% analyze modules first
[~,i,j]=intersect(meta.samples_for_analysis.Samplename,meta.pathways.sample);


% construct normalized matrix of pathway abundances
samples = meta.samples_for_analysis(i,:);
pathway_matrix = meta.pathways.P(:,j);
pathway_matrix = cell2mat(cellfun(@(x) x./sum(x), num2cell(pathway_matrix,1),'uni',0));

% for each pathway, compute the pearson correlation coefficient
for i = 1:length(meta.pathways.id)
    [r(i),p(i)]=corr(samples.Actualdepth,pathway_matrix(i,:)');
end

% find significant pathways using an FDR q-value of 0.05 (Benjamini &
% Hochberg procedure)
[h,~,~,padj]=fdr_bh(p,0.05);
pathways = meta.pathways.id(h);
[r_sorted,i]= sort(r(h),'descend');


% make a figure showing the pathways and their correlation with depth
barh(r_sorted);
set(gca,'YTick',1:length(r_sorted));
set(gca,'YTickLabel',pathways(i));
set(gca,'Fontsize',6.8);
ylabel('significant pathways')
xlabel('pearsons correlation coefficient')

% 
% g = table()
% g.pathways = meta.pathways.id'
% g.pearson_correlation_coeff = r'
% g.pearson_p_val = p'
% g.adjusted_pval = padj'
% g.significant = h'

