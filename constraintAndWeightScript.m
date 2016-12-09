for i = 1:length(subA)
    [m,n] = size(subB(i).B);
    subQ(i).Q = 1e-8*eye(m+n);
    subW(i).W = 1e3*eye(length(subA(i).A)/2);
    subR(i).R = -eye(size(subW(i).W));
    subV(i).V = 5e-2*eye(m+n);
end
for i = 1:length(subA)
    [n,~] = size(subH(i).H);
    subO_c(i).O_c = 15*ones(n,1);
    subI_c(i).I_c = 50;
end

[m,n] = size(Bhat);
Aaug = [Ahat,Bhat;zeros(n,m),-eye(n)];
Baug = [zeros(size(B));eye(n)];
[~,r] = size(Dhat);
Daug = [Dhat;zeros(n,r)];