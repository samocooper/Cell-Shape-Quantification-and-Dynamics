%% 

% -----------------Load Data------------------
cell_sample = csvread('Cell_sample.csv');
cell_sample = cell_sample(:,1:15);
[sample_size,feature_number] = size(cell_sample);
options = statset('MaxIter',500);


%%
% -------- Fit GMMs to each feature with number decided in opencv -------
feature_components = [2 1 1 2 1 1 1 1 2 1 1 1 2 2 2]; % used openCV in C++ to trial gmm number to each feature due to matlab crashing repeatedly
x_val_num = 10;
x_val = sample_size/x_val_num;
reduced_sample_set = [];
final_dist = [];
for i = 1:feature_number
    means = zeros(feature_components(i),1);
    covs = zeros(1,1,feature_components(i));
    for j = 1:x_val_num
        reduced_sample = cell_sample;
        for k = 1:x_val
            x = 10*x_val-k-j;
            reduced_sample(x,:) = [];
        end
        dist = gmdistribution.fit(reduced_sample(:,i),feature_components(i),'Options',options,'Replicates',10);
        means = means + dist.mu;
        covs = covs + dist.Sigma;
    end
    means = means/x_val_num;
    covs = covs/x_val_num;
    final_dist{i} = gmdistribution(means,covs);
    clustered_feature = cluster(final_dist{i},cell_sample(:,i));
    clustered_data(:,i) = clustered_feature-1;
end

%%
% ------- Determine cluster number with k means -------
cluster_qual = zeros(1,10);
x_val_num = 50;
for i = 1:10
    i
    for j = 1:x_val_num
    clustered_data_copy = clustered_data;
        for k = 1:x_val
            x = 10*x_val-k-j;
            clustered_data_copy(x,:) = [];
        end 
        [k_clusters,c] = kmeans(clustered_data_copy,i,'emptyaction','singleton','dist','Hamming','Replicates',50,'options',options);
        [sil,h1] = silhouette(clustered_data_copy,k_clusters,'Hamming');
        cluster_qual(i) = cluster_qual(i) + mean(sil);
    end
end
cluster_qual = cluster_qual/10;

%%
% ensure value drop off at high cluster number i.e. 10
temp = cluster_qual;
while temp(10)-min(temp) > 0
    diff = (temp(10)-min(temp))/10;
    temp = temp - (1:11)*(diff+0.001*temp(10));      
end
plot(temp)

%%
clusters = 5;
[index,centers] = kmeans(clustered_data,clusters,'emptyaction','singleton','dist','Hamming','Replicates',1000,'options',options);

%%
% ------------------Find files-------------------

names = dir('tracked_cells');
fnum = length(names);
fnames = {};
fcount = 0;
for i = 1:fnum
    if strfind(names(i).name,'csv')>0
        fcount = fcount+1;
        file_path = strcat('tracked_cells/',names(i).name);
        fnames{fcount} = file_path;       
    end
end
%% 

% ------------------Run model on test conditions----------------

shape_dist = zeros(fcount,clusters);
norm_shape_dist = zeros(fcount,clusters);

for i = 1:fcount
    
    tc = csvread(fnames{i});
    tc = tc(:,1:15);
    tc_size = length(tc);
    clustered_tcs = zeros(tc_size,feature_number);
    
    for j = 1:feature_number
        clustered_feature = cluster(final_dist{j},tc(:,j));
        clustered_tcs(:,j) = clustered_feature -1;
    end
    [dist_tc,index_tc] = min(pdist2(clustered_tcs,centers),[],2);
    
    for j = 1:clusters
        
        shape_dist(i,j) = sum(index_tc == j);
        norm_shape_dist(i,j) = sum(index_tc == j)/tc_size;
        %color = rand(1,3);
        %parallelcoords(tc((index_tc == j),:),'Color',color)
    end
    
end
bar(norm_shape_dist,'stacked')
