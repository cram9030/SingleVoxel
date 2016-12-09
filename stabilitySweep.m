function [contStab,contEig,contStiff,Kgain,Aaug,Baug,Daug] = stabilitySweep(m,k,c,num,divRange,weightRange,Vweight,Wweight)

contStab = zeros(length(divRange),length(weightRange),length(weightRange));
contEig = zeros(length(divRange),length(weightRange),length(weightRange));
contStiff = zeros(length(divRange),length(weightRange),length(weightRange));
Kgain(1,1,1).K = [];
Aaug(1).A = [];
Baug(1).B = [];
Daug(1).D = [];

for i = 1:length(divRange)
    ControlLocs = num/divRange(i):num/divRange(i):num;
    SensorLocs = [ControlLocs,ControlLocs-1];
    SensorLocs = sort(SensorLocs);
    [A,B_dist,B,C] = singleSpringSystem(m,k,c,num,SensorLocs,ControlLocs);
    [Ahat,Dhat,Bhat,Hhat,Mhat] = diagonalizeA(A,B_dist,B,C,C);
    [~,subA,~,subB,~,subD,~,subM,~,subH] = expandA(Ahat,Bhat,Dhat,Mhat,Hhat,ControlLocs,num);
    [~,n] = size(Bhat);
    [~,m] = size(Ahat);
    [~,r] = size(Dhat);
    Aaug(i).A = [Ahat,Bhat;zeros(n,m),-eye(n)];
    Baug(i).B = [zeros(size(B));eye(n)];
    Daug(i).D = [Dhat;zeros(n,r)];
    for k = 1:length(subA)
        subW(k).W = Wweight*eye(length(subA(k).A)/2);
        subI_c(k).I_c = 50;
        [n,~] = size(subH(k).H);
        subO_c(k).O_c = .15*ones(n,1);
        [m,n] = size(subB(k).B);
        subV(k).V = Vweight*eye(m+n);
    end
    for j = 1:length(weightRange)
        for k = 1:length(weightRange)
            for l = 1:length(subA)
                subR(l).R = weightRange(j)*eye(size(subW(l).W));
                [m,n] = size(subB(l).B);
                subQ(l).Q = weightRange(k)*eye(m+n);
            end
            [~,~,~,K,~] = DecentralizedLMI(subA,subD,subB,subH,subM,subW,subV,subQ,subR,subI_c,subO_c);
            Kgain(i,j).K = K;
            contEig(i,j) = max(real(eig(Aaug(i).A+Baug(i).B*K)));
            contStiff(i,j) = min(real(eig(Aaug(i).A+Baug(i).B*K)));
            if max(real(eig(Aaug(i).A+Baug(i).B*K)))<0
                contStab(i,j) = 1;
            else
                contStab(i,j) = 0;
            end
        end
    end
end