%Compute random component
function [M1,out,thrsh,random_component] = compute_M1_w_RMT(data);
%% Step 1A: Compute random expectation 

% Read package information for CP_ALS here
% (https://www.tensortoolbox.org/cp_als_doc.html#2).
% For determining the rank, we used the lowest value in a dimemsion within
% your tensor (e.g., a tensor of 3 people by 4 CAZymes by 10 timepoints
% (X=[3x4x10]) use 3 as the input rank).
prompt = 'In computing random expectation, how many iterations desired?';
iter = input(prompt);
clear prompt

prompt = 'In computing random expectation, what is desired rank?';
rank = input(prompt);
clear prompt

prompt = 'Do you want to compute covariance matrix?';
cov_decision = input(prompt); %1 is yes, 2 is no
clear prompt

% Start with 10, 20, 50, 100 iterations, and check histograms to determine distribution
% of the data (you can decide which to use for reporting the data).
prompt = 'In computing distribution of projection approximation, how many iterations of the bootstrap are desired?';
boots_dat = input(prompt);
clear prompt

sigma_out = zeros(rank,iter);
parfor k=1:iter; %Loop through number of iterations
    randomized_mat = zeros(length(data(:,1,1)),length(data(1,:,1)),length(data(1,1,:))); %Initialize randomized matrix
    M1_rand = [];
    
    if cov_decision == 1;
        for i = 1:length(data(1,1,:)); %Loop over 3rd dimension
            a = find(sum(data(:,:,i)') > 0); %Isolate rows that are sampled
            tmp = zeros(a,length(data(1,:,1)));
            tmp = data(a,i,:);
            tmp2 = zeros(length(tmp(:,1)),length(tmp(1,:)));
            
            for j=1:length(tmp(:,1));
                [m,n] = size(tmp);
                tmp2(j,randperm(n)) = tmp(j,:);
            end;
            
            randomized_mat(:,:,i) = cov(tmp2);
        end;
    else
        for i=1:length(data(1,1,:));
            tmp = zeros(length(data(:,1,1)),length(data(1,:,1)));
            tmp = data(:,:,i);
            tmp2 = zeros(length(tmp(:,1)),length(tmp(1,:)));
            
            for j=1:length(tmp(:,1));
                [m,n] = size(tmp);
                tmp2(j,randperm(n)) = tmp(j,:);
            end;
            
            randomized_mat(:,:,i) = tmp2;
        end;
    end;
    
    randomized_mat = tensor(randomized_mat);
    M1_rand = cp_als(randomized_mat,rank);
    sigma = M1_rand.lambda(1:rank);
    sigma = sigma./sum(sigma); %convert to variance
    sigma_out(:,k) = sigma;
end;

avg_sigma = mean(sigma_out'); %Average core tensor singular values variance for N iterations of randomized tensor formulation
figure; histogram(sigma_out,50); %Plots distribution of core tensor singular values
title('Distribution of random singular values');

X = tensor(data);
thrsh = avg_sigma(1); %First random averaged tensor component
clear i
clear sigma
counter = 0;

parfor l=1:boots_dat; %Loop through number of bootstrap iterations
    M1_boot_mean = [];
    M1_boot_med = []
    M1_boot_SD = [];
       for i=1:rank;
        counter = counter + 1;
        M1 = cp_als(X,i);
        sigma = M1.lambda;
        sigma = sigma./sum(sigma); %convert to variance
        out(i) = min(sigma);    
        end;
        M1_boot_mean = mean(M1(;
end;



figure; plot(out); title('Tensor singular value variance vs rank specified');

random_component = find(out > thrsh);
random_component = random_component(end);
clear M1
M1 = cp_als(X,random_component);

            
        
    
    
 