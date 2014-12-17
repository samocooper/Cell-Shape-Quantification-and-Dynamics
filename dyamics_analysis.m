% ---- Find files define arrays matrices etc. ----
orig_cell_sample = csvread('Cell_sample.csv');
warning('off','stats:kmeans:EmptyClusterRep')
names = dir('tracked_cells');

fnum = length(names);
fnames = {};
fcount = 0;
prev_gene = [];
gene_index = [];
gene_count = 0;
gene_names = {};

% ---- Link gene knockdowns to matrices ----

for i = 1:fnum
    if strfind(names(i).name,'csv')>0        
        fcount = fcount+1;
        file_path = strcat('tracked_cells/',names(i).name);
        fnames{fcount} = file_path;
        a = strfind(names(i).name,'_');
        gene = names(i).name(1:a-1);
        gene_names{fcount} = gene;
        if strcmp(gene,prev_gene)
            gene_index(fcount)= gene_count;            
        else
            prev_gene = gene;
            gene_count = gene_count + 1;            
            gene_index(fcount)= gene_count;
        end
                    
    end
end

%%

% ---- Open files from previous step ----

fileID1 = fopen('feature_set.txt','r');
feature_set = textscan(fileID1,'%d');
feature_set = cell2mat(feature_set);

fileID2 = fopen('centers.txt','r');
centers = textscan(fileID2,'%f %f %f %f %f');
centers = cell2mat(centers)';

fileID3 = fopen('means.txt','r');
means = textscan(fileID3,'%f');
means = cell2mat(means)';

dist = csvread('dist.csv');

%% 

clusters = 5
shape_change = zeros(clusters);
temp_means = means(feature_set);
tc_dynamics = zeros(fcount,clusters^2);

Mat_dists = zeros(clusters,clusters,fcount);

for h = 1:fcount
    temp_shape_change = zeros(clusters);
    fid = fopen(fnames{h});
    tcells = {};
    fline = fgetl(fid);
    tcell = str2num(fline(1:end-1));
    tcell = tcell(feature_set);
    i = 0;
    while feof(fid) == false
        fline = fgetl(fid);
        if isempty(fline)
             i = i+ 1;
             tcells{i} = tcell;
             tcell = [];
        else
            temp_tcell = str2num(fline(1:end-1));
            temp_tcell = temp_tcell(feature_set);
            tcell = [tcell; temp_tcell];
        end
    end
    cell_count = length(tcells);
    for i = 1:cell_count
        [time_points,f] = size(tcells{i});
        tc_mean_mat = repmat(temp_means,time_points,1);
        
        % ---- Convert cells to binary values ----
        
        clustered_cell = tcells{i} > tc_mean_mat;
        clustered_cell = +clustered_cell;
        
        % ---- Asign to a shape cluster ----
        
        [dist_tc,index_tc] = min(pdist2(clustered_cell,centers),[],2);
        
        % ---- Has shaped changed since the last timepoint ----        
        
        for j = 1:time_points-1
            temp_shape_change(index_tc(j),index_tc(j+1)) = temp_shape_change(index_tc(j),index_tc(j+1)) + 1;
        end

    end
       
    temp_shape_dist_mat = repmat(dist(h,:),clusters,1)';
    norm_shape_change = temp_shape_change./temp_shape_dist_mat;
    
    Mat_dists(:,:,h) = norm_shape_change;
end
%%
% ---- Avergae shape transition matrices over repeats ----
Avg_mat_dists = zeros(clusters,clusters,19);
Temp_mat = Mat_dists(:,:,1);
unique_gene_names = cell(1,19);
reps = 1;
count = 1;
for i = 2:fcount
    if strcmp(gene_names{i},gene_names{i-1})
        Temp_mat = Temp_mat + Mat_dists(:,:,i);
        reps = reps + 1;
    else
        Avg_mat_dists(:,:,count) = Temp_mat/reps;
        reps = 1;
        unique_gene_names{count} = gene_names{i-1};
        count = count + 1;
        Temp_mat = Mat_dists(:,:,i);
    end
end
%%
ident = zeros(size(Avg_mat_dists));
for i = 1:19
    ident(:,:,i) = eye(5,5);
end

diag = ident.*Avg_mat_dists;
off_diag = (ones(size(Avg_mat_dists))-ident).*Avg_mat_dists;
diag = sum(sum(diag),2);
off_diag = sum(sum(off_diag),2);
dynamic_measure = off_diag(:)./diag(:);
[dynamic_measure_sorted, Ind] = sort(dynamic_measure);

%%
% ---- Plot shape transition matrices ----
for i = 1:19
    subplot(4,5,i)
    imagesc(log(Avg_mat_dists(:,:,Ind(i))))
    colormap(cool)
    temp_val = num2str(dynamic_measure(Ind(i)),3);
    ti = strcat(unique_gene_names{Ind(i)},{' '},temp_val);
    title(ti)
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    %shape_change = shape_change + norm_shape_change;
    %tc_dynamics(h,:) = reshape(norm_shape_change,1,[]);    
end

%%

% ---- Cluster shape distributions (this is performed here as opposed to
% other matlab file as data is loaded from csv/txt files) ----

dist_avg = zeros(19,5);
Temp_dist_mat = dist(1,:);
count = 1;
reps = 1;
for i = 2:fcount
    if strcmp(gene_names{i},gene_names{i-1})
        Temp_dist_mat = Temp_dist_mat + dist(i,:);
        reps = reps + 1;
    else
        dist_avg(count,:) = Temp_dist_mat/reps;
        reps = 1;
        count = count + 1;
        Temp_dist_mat = dist(1,:);
    end
end

dist_avg = (dist_avg-0.25)*2;
clustergram(dist_avg,'RowLabels', unique_gene_names,'Colormap', 'redbluecmap');