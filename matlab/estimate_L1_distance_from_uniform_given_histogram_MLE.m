function d = estimate_L1_distance_from_uniform_given_histogram_MLE(hist_vec, S)
% estimates L1 distance from uniform distribution given histogram and support set size S
% uses plug-in (MLE) estimator
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
%     >> d = estimate_L1_distance_from_uniform_given_histogram_MLE([9,3,2,1,1], 8)

Z = sum(hist_vec);
d = sum(abs((hist_vec./Z)-1/S));

