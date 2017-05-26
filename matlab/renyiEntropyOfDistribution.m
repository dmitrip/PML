function H = renyiEntropyOfDistribution(vec,alpha,base)
% computes Renyi entropy of vector vec, parameter alpha, default base 2

if nargin == 2,
    base = 2;
end

vec = vec./sum(vec);

H = (log(sum(vec.^alpha))/(1-alpha)) / log(base);