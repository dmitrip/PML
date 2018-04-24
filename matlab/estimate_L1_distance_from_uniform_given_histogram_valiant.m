function d = estimate_L1_distance_from_uniform_given_histogram_valiant(hist_vec, S)
% estimates L1 distance from uniform distribution given histogram and support set size S
% uses Valiant, Valiant: http://theory.stanford.edu/~valiant/papers/VV_stoc11.pdf
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
%     >> d = estimate_L1_distance_from_uniform_given_histogram_valiant([9,3,2,1,1], 8)

%h1 = hist(samp,double(min(samp)):double(max(samp)));
h1 = hist_vec;
fingerprint=hist(h1,0:max(h1));
fingerprint = fingerprint(2:end);
[h,x]=unseen(fingerprint);
%valiant = -h*(x.*log(x))' / log(2);
%K_est = sum(h);
d = h*abs(x - 1/S)';


function [histx,x] = unseen(f)
 % Input: fingerprint f, where f(i) represents number of elements that
 % appear i times in a sample. Thus sum_i i*f(i) = sample size.
 % File makeFinger.m transforms a sample into the associated fingerprint.
 % Output: approximation of 'histogram' of true distribution. Specifically,
 % histx(i) represents the number of domain elements that occur with
 % probability x(i). Thus sum_i x(i)*histx(i) = 1, as distributions have
 % total probability mass 1.
 %
 % An approximation of the entropy of the true distribution can be computed
 % as: Entropy = (-1)*sum(histx.*x.*log(x))

 f=f(:)';
 k=f*(1:size(f,2))'; %total sample size


 %%%%%%% algorithm parameters %%%%%%%%%%%
 gridFactor = 1.1; % the grid of probabilities will be geometric, with this ratio.
 % setting this smaller may slightly increase accuracy, at the cost of speed
 alpha = .5; %the allowable discrepancy between the returned solution and the "best" (overfit).
 % 0.5 worked well in all examples we tried, though the results were nearly indistinguishable
 % for any alpha between 0.25 and 1. Decreasing alpha increases the chances of overfitting.
 xLPmin = 1/(k*max(10,k)); % minimum allowable probability.
 % a more aggressive bound like 1/k?.5 would make the LP slightly faster,
 % though at the cost of accuracy
 maxLPIters = 1000; % the 'MaxIter' parameter for Matlab's 'linprog' LP solver.
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 % Split the fingerprint into the 'dense' portion for which we
 % solve an LP to yield the corresponding histogram, and 'sparse'
 % portion for which we simply use the empirical histogram
 x=0;
 histx = 0;
 fLP = zeros(1,max(size(f)));
 for i=1:max(size(f))
 if f(i)>0
 wind = [max(1,i-ceil(sqrt(i))),min(i+ceil(sqrt(i)),max(size(f)))];
 if sum(f(wind(1):wind(2)))<2*sqrt(i)
 x=[x, i/k];
 histx=[histx,f(i)];
 fLP(i)=0;
 else
 fLP(i)=f(i);
  end
 end
 end

 % If no LP portion, return the empirical histogram
 fmax = max(find(fLP>0));
 if min(size(fmax))==0
 x=x(2:end);
 histx=histx(2:end);
 return;
 end

 % Set up the first LP
 LPmass = 1 - x*histx'; %amount of probability mass in the LP region

 fLP=[fLP(1:fmax), zeros(1,ceil(sqrt(fmax)))];
 szLPf=max(size(fLP));

 xLPmax = fmax/k;
 xLP=xLPmin*gridFactor.^(0:ceil(log(xLPmax/xLPmin)/log(gridFactor)));
 szLPx=max(size(xLP));

 objf=zeros(szLPx+2*szLPf,1);
 objf(szLPx+1:2:end)=1./(sqrt(fLP+1)); % discrepancy in ith fingerprint expectation
 objf(szLPx+2:2:end)=1./(sqrt(fLP+1)); % weighted by 1/sqrt(f(i) + 1)

 A = zeros(2*szLPf,szLPx+2*szLPf);
 b=zeros(2*szLPf,1);
 for i=1:szLPf
 A(2*i-1,1:szLPx)=poisspdf(i,k*xLP);
 A(2*i,1:szLPx)=(-1)*A(2*i-1,1:szLPx);
 A(2*i-1,szLPx+2*i-1)=-1;
 A(2*i,szLPx+2*i)=-1;
 b(2*i-1)=fLP(i);
 b(2*i)=-fLP(i);
 end

 Aeq = zeros(1,szLPx+2*szLPf);
 Aeq(1:szLPx)=xLP;
 beq = LPmass;


 options = optimset('MaxIter', maxLPIters,'Display','off');
 for i=1:szLPx
 A(:,i)=A(:,i)/xLP(i); %rescaling for better conditioning
 Aeq(i)=Aeq(i)/xLP(i);
 end
 [sol, fval, exitflag, output] = linprog(objf, A, b, Aeq, beq, ...
zeros(szLPx+2*szLPf,1), Inf*ones(szLPx+2*szLPf,1),[], options);
 if exitflag==0
    'maximum number of iterations reached--try increasing maxLPIters'
 end
 if exitflag<0
    'LP1 solution was not found, still solving LP2 anyway...'
 exitflag
 end

 % Solve the 2nd LP, which minimizes support size subject to incurring at most
 % alpha worse objective function value (of the objective function in the
 % previous LP).
 objf2=0*objf;
 objf2(1:szLPx) = 1;
 A2=[A;objf']; % ensure at most alpha worse obj value
 b2=[b; fval+alpha]; % than solution of previous LP
 for i=1:szLPx
 objf2(i)=objf2(i)/xLP(i); %rescaling for better conditioning
 end
 [sol2, fval2, exitflag2, output] = linprog(objf2, A2, b2, Aeq, ...
beq, zeros(szLPx+2*szLPf,1), Inf*ones(szLPx+2*szLPf,1),[], ...
options);

 if not(exitflag2==1)
 'LP2 solution was not found'
 exitflag2
 end


 %append LP solution to empirical portion of histogram
 sol2(1:szLPx)=sol2(1:szLPx)./xLP'; %removing the scaling
 x=[x,xLP];
 histx=[histx,sol2'];
 [x,ind]=sort(x);
 histx=histx(ind);
 ind = find(histx>0);
 x=x(ind);
 histx=histx(ind);
 
 
 



     
    
