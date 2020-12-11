function [ probCop,probEmp,probInd ] = copulaPredictionDistr( copulaFamily , th , binr1 , binr2 , nSamples)
% Function for predicting the joint probability disctribution between two
% neurons, given their marginals and using the copula model
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
% OUTPUTS:
%
% probCop = cell array over the timebins with the copula predicted joint distribution 
%
% probEmp = cell array over the timebins with the empirical joint distribution 
%
% probInd = cell array over the timebins with the independent joint distribution 
%

nMax = max( [binr1(:);binr2(:)]);
[R T] = size(binr1);

%predCov = zeros([T 1]);

probCop = cell([T 1]);
[probCop{:}] = deal(zeros(nMax+1));

probEmp = cell([T 1]);
[probEmp{:}] = deal(zeros(nMax+1));

probInd = cell([T 1]);
[probInd{:}] = deal(zeros(nMax+1));

samples = copularnd(copulaFamily,th,nSamples);

for tt=1:T
    tt;
    p1 = histcounts(binr1(:,tt),0:1:(nMax+1));
    p1=p1/sum(p1);
    p2 = histcounts(binr2(:,tt),0:1:(nMax+1));
    p2=p2/sum(p2);
    
    %% Independent
    probInd{tt} = p2 .* p1';

    %% Empirical
    probEmp{tt} = histcounts2(binr1(:,tt),binr2(:,tt),0:1:(nMax+1),0:1:(nMax+1));
    probEmp{tt} = probEmp{tt}./sum( probEmp{tt}(:) );

    
    %% Copula
    c1 = cumsum(p1)';
    c2 = cumsum(p2)';
    
    spikeCountPair = zeros([nSamples 2]);
    
    for sc = nMax+1:-1:1
        spikeCountPair( samples( : , 1) <c1(sc)   ,1) = sc-1;
        spikeCountPair( samples( : , 2) <c2(sc)   ,2) = sc-1;
        
    end
    
    
    for n1 = 0:nMax
        for n2 = 0:nMax
            probCop{tt}(n1+1,n2+1) = sum( (spikeCountPair(:,1)==n1) & (spikeCountPair(:,2)==n2));
            
        end
    end
    
    probCop{tt} = probCop{tt}/nSamples;
    
    probCop{tt} = max( probCop{tt} ,1/nSamples);
end

end

