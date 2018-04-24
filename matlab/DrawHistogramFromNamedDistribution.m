function [hist_vec, p] = DrawHistogramFromNamedDistribution(distr_name,n,S,param)
% draws histogram from distribution
%
% defines distributions from Valiant and Valiant paper ("VV") - 
% "Estimating the unseen: an n/log(n)-sample estimator 
% for entropy and support size, shown optimal via new CLTs"
% by Gregory Valiant and Paul Valiant, STOC 2011
%
% defines distributions from JVHW estimator paper - 
% "Minimax Estimation of Functionals of Discrete Distributions"
% by Jiantao Jiao, Kartik Venkat, Yanjun Han, and Tsachy Weissman
% S parameter of distribution, support size
% 
% created by Dmitri S. Pavlichin on May 14, 2017
%
% Matlab version: R2015a
%
% Args:
%     * distr_name - string, distribution name
%     * n - integer, sample size
%     * S - support set size, has different meanings for different
%     distributions
%
% Returns:
%     * hist_vec - S by 1, empirical histogram
%     Optional:
%     * H - base 2 entropy of distribution (true, not estimated) from which
%     sampled
%
% Example:
%     to draw the histogram:
%         >> hist_vec = DrawFromDistrRef('UnifVV',1000,500);
%     to draw the histogram and compute the entropy:
%         >> [hist_vec, H] = DrawFromDistrRef('UnifVV',1000,500);
p = [];
switch distr_name
    case 'UnifJVHW' % same as UnifVV
        p = ones(S,1);
        hist_vec = DrawFromMultinomial(p,n);
    case 'ZipfJVHW' % same as ZipfVV
        p = 1./(1:S)';
        hist_vec = DrawFromMultinomial(p,n);
    case 'UnifVV' % same as UnifJVHW
        p = ones(S,1);
        hist_vec = DrawFromMultinomial(p,n);
    case 'MixUnifVV'
        p = [(5/(2*S)).*ones(S/5,1) ; (5/(8*S)).*ones(4*S/5,1)];
        hist_vec = DrawFromMultinomial(p,n);
    case 'MixUnifVV_v2'
        p = [(5/(8*S)).*ones(4*S/5,1) ; (5/(2*S)).*ones(S/5,1)];
        hist_vec = DrawFromMultinomial(p,n);
    case 'ZipfVV' % same as ZipfJVHW
        p = 1./(1:ceil(S))';
        hist_vec = DrawFromMultinomial(p,n);
    case 'Zipf2VV'
        p = 1./((1:ceil(S))'.^(0.6));
        hist_vec = DrawFromMultinomial(p,n);
    case 'Zipf1WuYang' % same as ZipfVV, ZipfJVHW
        p = 1./(1:S)';
        hist_vec = DrawFromMultinomial(p,n);
    case 'Zipf2WuYang'
        p = 1./((1:ceil(S))'.^(0.5));
        hist_vec = DrawFromMultinomial(p,n);
    case 'MixtureWuYang'
        p = 1./((1:ceil(S))'.^(0.5));
        p((S/2+1):S) = (1-2/S).^((1:(S/2))-1);
        % normalize so both halves have mass 1/2
        Z1 = sum(p(1:round(S/2)));
        Z2 = sum(p((round(S/2)+1):end));
        p(1:round(S/2)) = 0.5.*p(1:round(S/2))./Z1;
        p((round(S/2)+1):end) = 0.5*p((round(S/2)+1):end)./Z2;        
        hist_vec = DrawFromMultinomial(p,n);
    case 'GeomVV'
        hist_vec = random('geom',1/ceil(S),[n 1])+1; % Geom(n/2) in VV notation
        hist_vec = histc(hist_vec,1:max(hist_vec));
    case 'MixGeomZipfVV'
        % choose which samples from Zipf, Geom
        n1 = random('bino',n,0.5);
        
        % draw from Zipf(n/2) distribution in VV notation
        p_Zipf = 1./(1:ceil(S/2))'; % Zipf(n/2) in VV notation
        p_Zipf = p_Zipf./sum(p_Zipf);
        % draw sample
        hist_Zipf = DrawFromMultinomial(p_Zipf,n1);
        
        % draw from geom distribution
        hist_Geom = random('geom',1/ceil(S/2),[n - n1,1])+1; % Geom(n/2) in VV notation
        hist_Geom = histc(hist_Geom,1:max(max(hist_Geom),length(hist_Zipf)));
        
        % combine samples
        hist_Zipf = [hist_Zipf(:) ; zeros(length(hist_Geom)-length(hist_Zipf),1)];
        hist_vec = hist_Zipf + hist_Geom;
        
        % get entropy
        p = pdf('geo',(1:100*ceil(S/2))-1,1/ceil( S/2));
        p = p./sum(p);
        p_Zipf = [p_Zipf(:) ; zeros(length(p)-length(p_Zipf),1)];
        
        p = 0.5.*(p(:) + p_Zipf(:));
        
    case 'MixGeomZipfWuYang'
        p1 = 1./(1:(S/2));
        p1 = p1./sum(p1);
        p2 = (1-2/S).^((1:(S/2))-1);
        p2 = p2./sum(p2);
        p = [0.5.*p1(:) ; 0.5.*p2(:)];
        hist_vec = DrawFromMultinomial(p,n);
    otherwise
        warning(['unrecognized distribution name: ' distr_name]);
end

p = p./sum(p);
