function K_est = estimate_support_from_histogram_PML_approximate(hist_vec)
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
%     * hist_vec - (vector) histogram, vector of nonnegative integers
%
% Returns:
%     * K_est - (integer) estimate for support set size
%
% Example:

% estimate support set size
[~,F0_est,~,not_have_valid_K_est] = PMLdistributionApproximate(hist_vec);
K_est = sum(hist_vec > 0) + F0_est;

if not_have_valid_K_est
    warning('do not have valid estimate for support set size')
end

