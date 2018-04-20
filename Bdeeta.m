function B=Bdeeta(eta,nm,beta,flag)
neta=size(eta,1);
nb=2*nm;
B=zeros(neta,nb);
if (flag==0)
    B(:,1)=ones(neta,1);
    for i=1:nm-1
        alphai=exp(-i*beta)-1;
        B(:,i+1)=(1+alphai)*eta./(1+alphai*eta);
    end
    for i=1:nm
        alphai=exp(-(i-1)*beta)-1;
        B(:,nm+i)=eta./(1+alphai-alphai*eta);
    end
else
    B(:,1)=zeros(neta,1);
    for i=1:nm-1
        alphai=exp(-i*beta)-1;
        B(:,i+1)=(1+alphai)./((1+alphai*eta).^2);
    end
    for i=1:nm
        alphai=exp(-(i-1)*beta)-1;
        B(:,nm+i)=(1+alphai)./((1+alphai-alphai*eta).^2);
    end
end
