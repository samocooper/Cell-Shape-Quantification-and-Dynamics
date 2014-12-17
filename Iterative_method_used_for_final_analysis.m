% ---- Load pooled sample, find files for each knockdown ----
cell_sample = csvread('Cell_sample.csv');
warning('off','all')
names = dir('tracked_cells');

% ---- Create index linking matrices to gene knockdowns ----

fnum = length(names);
fnames = {};
fcount = 0;
prev_gene = [];
gene_index = [];
gene_count = 0;
gene_names = {};

for i = 1:fnum
    if strfind(names(i).name,'csv')>0 
        fcount = fcount+1;
        file_path = strcat('tracked_cells/',names(i).name);
        fnames{fcount} = file_path;
        a = strfind(names(i).name,'_');
        gene = names(i).name(1:a-1);
        if strcmp(gene,prev_gene)
            gene_index(fcount)= gene_count;            
        else            
            gene_count = gene_count + 1;             
            gene_names{gene_count} = gene;
            prev_gene = gene;                    
            gene_index(fcount)= gene_count;
        end
        
    end
end
%%

% ---- Allocate space for matrices that will be called in Main loop ----

feat_red_profileD = zeros(1,13);
feat_red_profileM = zeros(1,13);

dim = 15;
feat_set = (1:dim);
clusters = 5;

[cell_num,~] = size(cell_sample);

% ---- Turn Pooled sample into binary profile ----

col_mean = sum(cell_sample)/cell_num;
mean_mat = repmat(col_mean,cell_num,1);

binary_profile = cell_sample > mean_mat;
binary_profile = +binary_profile;

% ---- Turn each gene knocdown (TC) into binary profiles ----

tc_profiles = cell(fcount,1);
tc_sizes = zeros(fcount,1);
for i = 1:fcount
        tc = csvread(fnames{i});
        
        [tc_size,fnum] = size(tc);
        tc_mean_mat = repmat(col_mean,tc_size,1);
        
        binary_tc = tc > tc_mean_mat;
        binary_tc = +binary_tc;
        
        tc_profiles{i} = binary_tc;
        tc_sizes(i,1) = tc_size;
end
tc_sizes = repmat(tc_sizes,1,clusters);
iters = 13;

%%

% ---- Main loop ----

norm_dists = cell(iters,1);
feat_sets = cell(iters,1);
centers_sets = cell(iters,1);

for g = 1:iters % one feature is removed after each iteration
    
    tc_mean = zeros(1,dim);
    tc_dist = zeros(1,dim);
    temp_dists = cell(dim,1);
    temp_centers_sets = cell(iters,1);
    parfor h = 1:dim
        
        % ---- In this iteration remove one feature test out how it 
        % effects the distance between repeats compared to distance 
        % between gene knockdowns(TCs) ----
        
        % -- Determine cluster centers from pooled sample using k means --
        
        warning('off','all')
        
        feat_set_temp = feat_set;
        feat_set_temp(h) = [];
        
        binary_temp = binary_profile(:,feat_set_temp);
        
        [index,centers] = kmeans(binary_temp,clusters,'MaxIter',800,'emptyaction','singleton','dist','hamming','Replicates',600);
        
        cluster_size = zeros(1,clusters);
        for i = 1:clusters
            cluster_size(i) = sum(index == i)/cell_num;
        end
        
        % -- Run model on test conditions --

        shape_dist = zeros(fcount,clusters);
        
        for i = 1:fcount
            binary_tc = tc_profiles{i}(:,feat_set_temp);
            [dist_tc,index_tc] = min(pdist2(binary_tc,centers),[],2);
            
            for j = 1:clusters
                shape_dist(i,j) = sum(index_tc == j);
            end
        end
        
        norm_dist = shape_dist./tc_sizes;
        
        temp_dists{h} = norm_dist;
        temp_centers_sets{h} = centers;
        
        gene_dev = zeros(1,clusters);
        gene_dif = zeros(1,clusters);
        
        for i = 1:gene_count
            
            tc_in_gene = sum(gene_index==i);
            
            if (tc_in_gene>1)
                
                gene_mean = sum(norm_dist(gene_index == i,:))/tc_in_gene;                
                gene_dif = gene_dif + abs(gene_mean-cluster_size);
                
                gene_mean = repmat(gene_mean,tc_in_gene,1);
                gene_dev = gene_dev + sum(abs(norm_dist(gene_index == i,:)-gene_mean));
            end
        end        
        tc_dist(h) = sum(gene_dev);   
        tc_mean(h) = sum(gene_dif);     
    end
    
    % ---- Save results from each iteration for later plotting ----
    [~,ind] = max(tc_mean-tc_dist);
    
    [val1,~] = min(tc_dist);    
    [val2,~] = max(tc_mean);
    
    feat_red_profileD(g) = val1;
    feat_red_profileM(g) = val2;
    
    
    feat_set(ind) = [];
    
    feat_sets{g} =  feat_set;    
    norm_dists{g} = temp_dists{ind};    
    centers_sets{g} = temp_centers_sets{ind};
    
    dim = dim - 1; 
    
end
%%

% ---- Save infomation for dynamic profiles ----

fileID1 = fopen('feature_set_2.txt','w');
fprintf(fileID1,'%d ',feat_sets{8});
fclose(fileID1);

fileID2 = fopen('centers_2.txt','w');
fprintf(fileID2,'%f %f %f %f %f\r\n',centers_sets{8});
fclose(fileID2);

csvwrite('dist_2.csv',norm_dists{8});

fileID5 = fopen('means_2.txt','w');
fprintf(fileID5,'%f\r\n',col_mean);
fclose(fileID5);
%%
% -------------------- Plotting Results ----------------------

% -- Sort rows by cluster size so same cluster is generally on bottom for each
% iteration of a feature removal (since unsupervised clusters otherwise switch
% randomly i.e. cluster 1 from iteration 1 is cluster 4 in iteration 4. --

for i = 1:iters
    norm_dist = norm_dists{i};
    norm_dist = norm_dist';
    norm_dist = sortrows(norm_dist,1);  
    norm_dist = flipud(norm_dist);
    norm_dist = norm_dist';    
    norm_dists{i} = norm_dist;
end
%%

% ---- Plot shape cluster distribution for all TCs for chosen feature
% reduction number ----

figure
bar(norm_dists{8},'stacked')
%%
% ---- Plot reduction in variation ----

hold on
plot(feat_red_profileD(1:12)/10) % variation between gene knockdowns
plot(feat_red_profileM(1:12)/10) % Variation within gene knockdowns
hold off
%%

% ---- Plot changes in cluster ditribution for a given knockdown ----

% -- linnked with previous figure for paper figure --

figure
hold on
for i = 1:iters
    norm_dist_temp = norm_dists{i};
    colour = [0 0 0.6];
    scatter([i i i i],norm_dist_temp(31:34,1),60,colour,'.')
end
for i = 1:iters
    norm_dist_temp = norm_dists{i};
    colour = [0.2 0.5 0.8];
    scatter([i i i i],norm_dist_temp(31:34,2),60,colour,'.')
end
for i = 1:iters
    norm_dist_temp = norm_dists{i};
    colour = [0.5 1 0.3];
    scatter([i i i i],norm_dist_temp(31:34,3),60,colour,'.')
end
for i = 1:iters
    norm_dist_temp = norm_dists{i};
    colour = [ 1 0.5 0];
    scatter([i i i i],norm_dist_temp(31:34,4),60,colour,'.')
end
for i = 1:iters
    norm_dist_temp = norm_dists{i};
    colour = [0.8 0 0];
    scatter([i i i i],norm_dist_temp(31:34,5),60,colour,'.')
end
hold off
%%

% ---- Plot change in distirbutions over 6 feature removals ----

figure
for h = 1:6

norm_dist = norm_dists{h*2};
gene_mean = zeros(gene_count,clusters);
for i = 1:gene_count
    tc_in_gene = sum(gene_index==i);
    if (tc_in_gene>1)        
        gene_mean(i,:) = sum(norm_dist(gene_index == i,:))/tc_in_gene;
    else
        gene_mean(i,:) = norm_dist(gene_index == i,:);
    end
end
subplot(2,3,h)
bar(gene_mean,'stacked')
axis([0 21 0 1])
end
