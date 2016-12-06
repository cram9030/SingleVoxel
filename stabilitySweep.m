function [contStab,obsrvStab] = stabilitySweep(m,k,c,num,divRange,weightRange)

contStab = zeros(length(divRange),length(weightRange),length(weightRange));
obsrvStab = zeros(length(divRange),length(weightRange),length(weightRange));

for i = 1:length(divRange)
    for j = 1:length(weightRange)
        ControlLocs = num/divRange(i):num/divRange(i):num;
        [A,B_dist,B,C] = singleSpringSystem(m,k,c,num,2*ControlLocs,ControlLocs);
        [Ahat,Dhat,Bhat,Hhat,Mhat] = diagonalizeA(A,B_dist,B,C,C);
        [Atilde,subA,Btilde,subB,Dtilde,subD,Mtilde,subM,Htilde,subH] = expandA(Ahat,Bhat,Dhat,ControlLocs,num);
        for k = 1:length(subA)
            subW(k).W = weightRange(j)*eye(length(subA(k).A)/2);
            subI_c(k).I_c = 50;
            [n,~] = size(subH(k).H);
            subO_c(k).O_c = .15*ones(n,1);
        end
        [~,n] = size(Bhat);
        [~,m] = size(Ahat);
        [s,~] = size(Mhat);
        [~,r] = size(Dhat);
        Aaug = [Ahat,Bhat;zeros(n,m),-eye(n)];
        Baug = [zeros(size(B));eye(n)];
        Maug = [Mhat,zeros(s,n)];
        for k = 1:length(weightRange)
            for l = 1:length(subA)
                subV(l).V = weightRange(k);
            end
            [subL,subK,subG,subZ,K,L] = DecentralizedLMI(subA,subD,subB,subH,subM,subW,subV,subI_c,subO_c);
            if max(real(eig(Aaug+L*Maug)))<0
                obsrvStab(i,j,k) = 1;
            else
                obsrvStab(i,j,k) = 0;
            end
            if max(real(eig(Aaug+Baug*K)))<0
                contStab(i,j,k) = 1;
            else
                contStab(i,j,k) = 0;
            end
        end
    end
end