
clear


% add paths
addpath('../lib/cbrewer/')
addpath('../lib/fdr_bh/')

load('../data/meta_full.mat');
load('../KEGG/modules_AB_noCyano_80percentRepresent.mat')
% construct the module matricies for only heterotrophs:


% normalize the metagenomes
X = cell2mat(cellfun(@(x) x./sum(x),num2cell(meta.hetero.metagenomes,1),'uni',0));

% for each module, find the sum of gene relative abundances, creaing a
% module matrix

M.X = [];
for i = 1:length(rules)
    [~,z]=intersect(meta.hetero.metagenome_genes,rules(i).ko);
    M.X = [M.X;sum(X(z,:),1)];
end
M.module = {rules.Module};

% get the sample idicies for different filter sizes 
idx11 = meta.samples_info.Filter == 11;
idx22 = meta.samples_info.Filter == 0.22;


% for each module and filter size, compute the pearson correlation between module
% and depth, the mean of each module across all depths, and the slope of
% the best-fit line between module fraction and depth
for i = 1:length(M.module)
    [r,p]=corr(M.X(i,idx22)',meta.samples_info.Depth_m(idx22));
    R22(i) = r;
    P22(i) = p;
    M22(i) = mean(M.X(i,idx22));
    B22{i} = glmfit(M.X(i,idx22)',meta.samples_info.Depth_m(idx22));
    
    
    
    [r,p]=corr(M.X(i,idx11)',meta.samples_info.Depth_m(idx11));
    R11(i) = r;
    P11(i) = p;
    M11(i) = mean(M.X(i,idx11));
    B11{i} = glmfit(M.X(i,idx11)',meta.samples_info.Depth_m(idx11));

    
end
    



nitrogen_modules = {'M00175','M00531','M00530','M00529','M00528'};
sulfur_modules = {'M00176','M00596','M00595'};
%phos_modules = {'
[~,z]=intersect(M.module,nitrogen_modules);
%[~,z]=intersect(M.module,sulfur_modules);

%name = rules(strcmp({rules.Module},M.module(z(1)))).Name;
name = arrayfun(@(y) rules(strcmp({rules.Module},M.module(y))).Name,z,'uni',0);

idx11 = meta.samples_info.Filter == 11;
idx22 = meta.samples_info.Filter == 0.22;

%% for each nitrogen module, compute the association between depth and fraction
cmap = [cbrewer('qual','Set1',9)];
figure();
feat = 'Depth_m';
numMod = length(z);
b =1;
for k = 1:length(z)
    subplot(numMod,2,b);
    scatter(M.X(z(k),idx22),meta.samples_info.(feat)(idx22),[],cmap(2,:),'filled','MarkerEdgeColor','k');
    [r,p]=corr(M.X(z(k),idx22)',meta.samples_info.(feat)(idx22));
    if p<0.05
        hold on
        text(max(M.X(z(k),idx22)).*0.9,400,num2str(r));
        hold off
    end
    title(name(k));
    P(k,1) = p;
    
    b = b+1;
    subplot(numMod,2,b);
    scatter(M.X(z(k),idx11),meta.samples_info.(feat)(idx11),[],cmap(1,:),'filled','MarkerEdgeColor','k');
    [r,p]=corr(M.X(z(k),idx11)',meta.samples_info.(feat)(idx11));
    if p<0.05
        hold on
        text(max(M.X(z(k),idx11)).*0.9,400,num2str(r));
        hold off
    end
    b = b+1;
    P(k,2) = p;
    %title(name(k));
end



