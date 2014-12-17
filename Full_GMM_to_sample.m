%% 
% ----------------- Preprocess data---------------------
cell_sample = csvread('Cell_sample.csv');
[sample_size,feature_number] = size(cell_sample);
options = statset('MaxIter',800);
cell_sample = cell_sample(:,1:15);
cell_sample = zscore(cell_sample);
[~,cell_sample,latent] = princomp(cell_sample);
variance = cumsum(latent)./sum(latent); % PCs 1 to 4 account for 90% of variation this is used to fit clusters to reduce dimensionality of problem
%%

% -----------------Generate GMM of sample------------------

BIC = zeros(1,10);
obj = cell(1,10);
for k = 1:10
    obj{k} = gmdistribution.fit(cell_sample(:,1:4),k,'Options',options,'Replicates',50);
    BIC(k)= obj{k}.BIC;
end
%%

% -------------------Plotting--------------------
[minBIC,numComponents] = min(BIC);
plot(BIC)
numComponents



