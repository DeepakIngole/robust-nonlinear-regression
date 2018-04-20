function wienerModel=robust_nlreg(trainData,options)

nm=options.nm; 
beta=options.beta;
eta=options.eta;
eps=options.eps;
sigma=options.sigma;
[nData,nz]=size(trainData(:,1:end-1));
nb=2*nm;
np=nb+nz+nData;

Pi=[-Bdeeta(eta,nm,beta,1)*Psel(1:nb,np);
    Bdeeta(trainData(:,end),nm,beta,0)*Psel(1:nb,np)-...
    trainData(:,1:end-1)*Psel(nb+1:nb+nz,np)-...
    Psel(nb+nz+1:np,np);
    -Bdeeta(trainData(:,end),nm,beta,0)*Psel(1:nb,np)+...
    trainData(:,1:end-1)*Psel(nb+1:nb+nz,np)-...
    Psel(nb+nz+1:np,np);
    -Psel(nb+nz+1:np,np)];
Upsilon=[-eps*ones(numel(eta),1);
         sigma*ones(2*nData,1);
         zeros(nData,1)];

weightFactor=ones(size(trainData,1),1);
Omega=ones(1,nData)*diag(weightFactor)*Psel(nb+nz+1:np,np);
sparsity_ratio_prev=0;
for k=1:options.iter_max
	p_opt=cplexlp(Omega,Pi,Upsilon);
    weightFactor=1./(options.delta+p_opt(nb+nz+1:end,1));
    Omega=ones(1,nData)*diag(weightFactor)*Psel(nb+nz+1:np,np); 
    sparsity_ratio=numel(find(p_opt(nb+nz+1:end,1)<=options.delta))/numel(p_opt);
    if sparsity_ratio<=sparsity_ratio_prev
        try
            p_opt=p_opt_prev;
            sparsity_ratio=sparsity_ratio_prev;
        end
        break;
    else
        sparsity_ratio_prev=sparsity_ratio;
        p_opt_prev=p_opt;
    end
end

wienerModel.mu=p_opt(1:nb,1);
wienerModel.L=p_opt(nb+1:nb+nz,1);
wienerModel.nm=nm; 
wienerModel.beta=beta;
wienerModel.eta=eta;

