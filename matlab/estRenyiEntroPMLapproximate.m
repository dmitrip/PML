function [RenyiH_est, K_est] = estRenyiEntroPMLapproximate(samp,alpha,K)
% estimates Renyi entropy with parameter alpha of distribution given sample samp
% optionally given assumed true support set size K
%
% approximates pattern maximum likelihood (PML) distribution p
% where p \approx \arg \max_p' P(unordered histogram)
% Approximation sums over all permutations that mix within blocks of
% constant p
%
% Uses approximate PML distribution as plug-in for Renyi entropy functional
%
% created by Dmitri S. Pavlichin on May 16, 2017
%
% Matlab version: R2015a
%
% Args:
%     * samp - (vector) sample, vector of positive integers
%     * alpha - (double) Renyi entropy parameter, alpha >= 0, alpha != 1
% Optional args:
%     * K - (integer) assumed support set size.  
%         If K is not provided, then we attempt to estimate the support set size
%
% Returns:
%     * H_est - (double) estimate for Renyi entropy, log base 2
%     * K_est - (integer) estimate for support set size
%
% Example:

% get empirical histogram
hist_vec = int_hist(samp(:));

if nargin == 3
    K_est = K;
    p_approx = PMLdistributionApproximate(hist_vec, K);
    not_have_valid_K_est = false;
else
    % estimate support set size
    [p_approx,F0_est,~,not_have_valid_K_est] = PMLdistributionApproximate(hist_vec);
    K_est = sum(hist_vec > 0) + F0_est;
end

if not_have_valid_K_est
    warning('do not have valid estimate for support set size')
end

% plug-in to Renyi entropy functional
RenyiH_est = renyiEntropyOfDistribution(p_approx, alpha, 2);
