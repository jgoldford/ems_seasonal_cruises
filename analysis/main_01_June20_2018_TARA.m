% analysis of TARA oceans data
% June-20-2018
% Joshua Goldford

% load the data
load('../TARA/tara_data_parsed.mat')

% analyze modules first
[~,i,j]=intersect(meta.samples_for_analysis.Samplename,meta.modules.samples);

samples = meta.samples_for_analysis(i,:);
module_matrix = meta.modules.M(:,j);
module_matrix = cell2mat(cellfun(@(x) x./sum(x), num2cell(module_matrix,1),'uni',0));


nitrogen_modules = {'M00175','M00531','M00530','M00529','M00528'};
sulfur_modules = {'M00176','M00596','M00595'};
[~,nit,nit_j]=intersect(meta.modules.id,nitrogen_modules);
[~,sulf,suf_j]=intersect(meta.modules.id,sulfur_modules);
sulfur_modules = sulfur_modules(suf_j);
nitrogen_modules = nitrogen_modules(nit_j);

% make figures for nitrogen and sulfur modules w/ depth
figure();
for i = 1:length(nit)
    subplot(length(nit),1,i);
    scatter(samples.Actualdepth,module_matrix(nit(i),:)','k');
    title(nitrogen_modules{i});
    xlabel('depth [m]')
    ylabel('metagenome fraction');
    [r,p]=corr(samples.Actualdepth,module_matrix(nit(i),:)','type','pearson');
    if p < 0.05
        text(800,0.8.*max(module_matrix(nit(i),:)),strcat('r^2 = ',num2str(r^2,2)));
    end
        
end


figure();
for i = 1:length(sulf)
    subplot(length(sulf),1,i);
    scatter(samples.Actualdepth,module_matrix(sulf(i),:)','k');
    title(sulfur_modules{i});
    xlabel('depth [m]')
    ylabel('metagenome fraction');
    [r,p]=corr(samples.Actualdepth,module_matrix(sulf(i),:)','type','pearson');
    if p < 0.05
        text(800,0.8.*max(module_matrix(sulf(i),:)),strcat('r^2 = ',num2str(r^2,2)));
    end
        
end


