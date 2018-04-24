%% get performance of estimators of L1 distance from uniform distribution, known support set size
addpath('../matlab')

%% parameters
S = 10000; % support set size
n_list = round(10.^(2:1:5)); % sample size
distribution_list = {'UnifVV','MixUnifVV','ZipfVV','Zipf2VV'};
num_trials = 10;

%% list of estimators sharing common invocation
L1_distance_from_uniform_estimators = {'MLE','valiant','PML_approx'};

num_est = length(L1_distance_from_uniform_estimators);

est_mat = zeros(num_trials, num_est, length(n_list),length(distribution_list));
L1_distance_from_uniform_list = zeros(num_trials, num_est, length(n_list),length(distribution_list));

%% get data
tic
for n_ix = 1:length(n_list)
    n = n_list(n_ix);
    disp(['n = ' num2str(n)])
    for distr_ix = 1:length(distribution_list)
        distribution_name = distribution_list{distr_ix};
        disp(['distribution: ' distribution_name])
            
        [~, p] = DrawHistogramFromNamedDistribution(distribution_name, n, S);
        
        L1_distance_from_uniform_true = sum(abs((p./sum(p)) - 1/S));
        
        for trial_num = 1:num_trials
            disp([n_ix distr_ix trial_num])
            [hist_, ~] = DrawHistogramFromNamedDistribution(distribution_name, n, S);
            
            % estimate L1 distance from uniform
            est_num = 0;
            for est = L1_distance_from_uniform_estimators
                est_num = est_num + 1;
                eval(['est_mat(trial_num,est_num,n_ix,distr_ix) = ' ...
                    'estimate_L1_distance_from_uniform_given_histogram_' est{1} '(hist_,S);']);
                L1_distance_from_uniform_list(trial_num,est_num,n_ix,distr_ix) = L1_distance_from_uniform_true;
            end
        end
    end
end
toc

%% save to .mat
save data_L1_distance_from_uniform_performance_test_v3