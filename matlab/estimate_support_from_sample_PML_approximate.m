function K_est = estimate_support_from_sample_PML_approximate(samp)
% estimates support set size of distribution given sample samp
%
% approximates pattern maximum likelihood (PML) distribution p
% where p \approx \arg \max_p' P(unordered histogram)
% Approximation sums over all permutations that mix within blocks of
% constant p
%
% Uses approximate PML distribution as plug-in for support functional
%
% created by Dmitri S. Pavlichin on June 8, 2017
%
% Matlab version: R2015a
%
% Args:
%     * samp - (vector) sample, vector of positive integers
%
% Returns:
%     * K_est - (integer) estimate for support set size
%
% Example:

% get empirical histogram
hist_vec = int_hist(samp(:));

K_est = estimate_support_from_histogram_PML_approximate(hist_vec);

