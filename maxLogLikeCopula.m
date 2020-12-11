function [ th  ] = maxLogLikeCopula(copulaFamily,bins,tol,bounds)
% Fitting discrete time-dependent copula model by maximum log-likelihood
%
% INPUTS
%
% copulaFamily = 'Gaussian', 'Clayton', 'Frank', or 'Gumbel'
%
% bins = cell array for the binned spikes of the two neurons where
% bins{t}(repeat , neur) = spike count of cell 'neur' during the repetition
% 'repeat' at time 't'
%
% tol = telorance for the minimization
%
% bounds = [lower_bond , upper_bound] for the copula parameter
%
% OUTPUT
%
% th = inferred copula parameter

T = numel(bins);
nMax = max(arrayfun(@(t) max(bins{t}(:)),1:T) ); 

probEmp = cell([T 1]);
cs1 = zeros([T nMax+1]);
cs2 = zeros([T nMax+1]);

for tt=1:T
    
    probEmp{tt} = histcounts2(bins{tt}(:,1),bins{tt}(:,2),0:1:(nMax+1),0:1:(nMax+1));
    probEmp{tt} = probEmp{tt}./sum( probEmp{tt}(:) );
    
    p1 = histcounts(bins{tt}(:,1),0:1:(nMax+1));
    p1=p1/sum(p1);
    p2 = histcounts(bins{tt}(:,2),0:1:(nMax+1));
    p2=p2/sum(p2);
    
    cs1(tt,:) = cumsum(p1);
    cs2(tt,:) = cumsum(p2);
    
    %pmfs{tt} = {p1,p2};
end


ops=optimset('TolX',tol,'maxiter',50,'display','off','LargeScale','off');
[th,vval] = fminbnd(@(x) SumTime(copulaFamily,x,probEmp,cs1,cs2,nMax),bounds(1),bounds(2),ops);


end

%% SumTime

function ellTime = SumTime(copulaFamily,x,probEmp,cs1,cs2,nMax)

T = numel(probEmp);
ellTime = - sum( arrayfun( @(tt)  logLikeCopula( copulaFamily, x , probEmp{tt},cs1(tt,:),cs2(tt,:),nMax  )   ,1:T) );


end


%% logLikeCopula

function ell =  logLikeCopula( copulaFamily, x , probEmp_t,cs1_t,cs2_t,nMax  ) 
[u1,u2] = meshgrid(cs2_t,cs1_t); % BE CAREFUL !!!

y = copulacdf(copulaFamily,[u1(:),u2(:)],x);

y = reshape(y,[nMax+1 nMax+1]);


logP = zeros(nMax+1);

logP(1,:) = [y(1,1),diff(y(1,:))];
logP(:,1) = [y(1,1);diff(y(:,1))];

for n1=1:nMax
    for n2=1:nMax
        logP(n1+1,n2+1) = y(n1+1,n2+1) - y(n1,n2+1) - y(n1+1,n2) + y(n1,n2);        
    end
end

logP = log(logP);

index = probEmp_t>0;
ell = sum( probEmp_t(index) .* logP(index) );




end

