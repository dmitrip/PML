function d = estimate_L1_distance_from_uniform_given_histogram_PML_approx(hist_vec, S)
% estimates L1 distance from uniform distribution given histogram and support set size S
%
% approximates pattern maximum likelihood (PML) distribution p
% where p \approx \arg \max_p' P(unordered histogram)
% Approximation sums over all permutations that mix within blocks of
% constant p
%
% created by Dmitri S. Pavlichin on August 16, 2017
%
% Matlab version: R2015a
%
% Args:
%     * hist_vec - (vector) histogram, vector of nonnegative integers
%     * S - (integer) support set size.  Must have S >= length(hist_vec)
%
% Returns:
%     * d - (float) estimate of L1 distance from uniform distribution with support set size S
%
% Example:
%     >> d = estimate_L1_distance_from_uniform_given_histogram_PML_approximate([9,3,2,1,1], 8)

% estimate PML distribution
p_approx = PMLdistributionApproximate(hist_vec, S);

d = sum(abs(p_approx - 1/S));

