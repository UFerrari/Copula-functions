function [ predCov ] = copulaPrediction( copulaFamily , th , binr1 , binr2 , nSamples)
% Function for predicting the noise covariances between two neurons 
% given their marginals and using the copula model
%
% INPUTS:
%
% copulaFamily = 'Gaussian', 'Clayton', 'Frank', or 'Gumbel'
%
% th = paramter of the copula model
%
% binr1,bir2 = arrays with the spike counts across time and repetitions
% size(binr1) =  [n_repetition, n_timebin]
%
% nSamples = numebr of random samples used for generating synthetic joint spike trains 
%
% OUTPUT:
%
% predCov = model predicted noise covariance for all timebins



[R, T] = size(binr1);

predCov = zeros([T 1]);

%nSamples = 50000;


nMax = max( [binr1(:);binr2(:)]);
for tt=1:T
    tt;
    p1 = histcounts(binr1(:,tt),0:1:(nMax+1));
    p1=p1/sum(p1);
    p2 = histcounts(binr2(:,tt),0:1:(nMax+1));
    p2=p2/sum(p2);
    
    c1 = cumsum(p1)';
    c2 = cumsum(p2)';
    
    spikeCountPair = zeros([nSamples 2]);
    samples = copularnd(copulaFamily,th,nSamples);
    for sc = nMax+1:-1:1
        spikeCountPair( samples( : , 1) <c1(sc)   ,1) = sc-1;
        spikeCountPair( samples( : , 2) <c2(sc)   ,2) = sc-1;
        
    end
    
    predCov(tt) = mean(prod(spikeCountPair,2)) - prod(mean(spikeCountPair)) ;
    
end

end

