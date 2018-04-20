function P=Psel(i,N)
P=zeros(length(i),N);
for k=1:length(i)
    P(k,i(k))=1;
end
