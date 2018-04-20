function [wienerModel,weightFactor]=wiener_train(wienerModel,trainData)
weightFactor=weight_initialize(trainData);
wienerModel=wiener_optimize(wienerModel,trainData,weightFactor); 
for i=1:3
    weightFactor=weight_tune(weightFactor,wienerModel);
    wienerModel=wiener_optimize(wienerModel,trainData,weightFactor);
end
wienerModel=rmfield(wienerModel,'err');

function weightFactor=weight_initialize(trainData)
l=size(trainData,1);
weightFactor=ones(l,1);
for k=1:l
    try
        stData=evalin('base','stData');
        for i=1:size(stData,1)
            if norm(trainData(k,:)-stData(i,:))<=1e-2 
                weightFactor(k)=2e1;
            end
        end 
    end
end

function weightFactor=weight_tune(weightFactor,wienerModel)
err=wienerModel.err;
c1=2.5;c2=3;
s=iqr(err)/2/.6745;
for k=1:length(err)
    if err(k)/s >=c1 && err(k)/s <=c2
        weightFactor(k)=(c2-err(k)/s)/(c2-c1);
    elseif err(k)/s >c2
        weightFactor(k)=1e-2*weightFactor(k);
    end
end

function wienerModel=wiener_optimize(wienerModel,trainData,weightFactor)

% Identification parameters
nm=wienerModel.nm; 
beta=wienerModel.beta;
eta=wienerModel.eta;
eps=wienerModel.eps;

%
[~,nz]=size(trainData(:,1:end-1));
nb=2*nm;
np=nb+nz;

% Linear least squares matrices
Omega=[Bdeeta(trainData(:,end),nm,beta,0) -trainData(:,1:end-1)];
neta=length(eta);
Pi=[-Bdeeta(eta,nm,beta,1)  zeros(neta,nz)];
Upsilon=-eps*ones(neta,1);

% QP matrices and preprocessing
H=Omega'*diag(weightFactor)*Omega;
H=.5*(H+H');
i=-30;
while 1
    [~,p]=chol(H+(10^i)*eye(np));
    if p==0
        i=i+5;
        H=H+(10^i)*eye(np);
        break;
    else
        i=i+1;
    end
end
f=zeros(1,np)';

% Solve QP
try
    p_opt=cplexqp(H,f,Pi,Upsilon);
catch
    p_opt=cplexqp(H+10^(-3)*eye(np),f,Pi,Upsilon);
end
mu=p_opt(1:nb);
L=p_opt(nb+1:nb+nz);

% Model
wienerModel.mu=mu;
wienerModel.L=L;

wienerModel.err=abs(Omega*p_opt);

