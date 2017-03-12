%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,G,Z] = LMIDesign(A,D,B,H,M,W,V,Q,R,I_c,O_c)

yalmip('clear');
[m,n] = size(B);
[~,r] = size(D);
[s,~] = size(M);
Aalpha = [A+7e0*eye(size(A)),B;zeros(n,m),zeros(n)];
A = [A,B;zeros(n,m),-eye(n)];
B = [zeros(size(B));eye(n)];
D = [D;zeros(n,r)];
M = [M,zeros(s,n);zeros(n,m),eye(n)];
C = [H,zeros(length(O_c),length(I_c));zeros(length(I_c),m),eye(length(I_c))];

Y = sdpvar(length(A),length(A));
Z = sdpvar(length(A),length(A));
G = sdpvar(n,m+n);

%estimator
%Y = care(A'+0.01*eye(size(A)),M',D*W*D'+0*eye(size(D*W*D')),V,zeros(size(M')),eye(size(A)));

%solve for observer gain
%Y = double(Y);
%L = Y*M'*inv(V);

%Controller
F = [[Z*Aalpha'+Aalpha*Z+G'*B'+B*G+Q  D*sqrt(W);sqrt(W)'*D' R]<0];
F = F + [Z-V> 0]; % %Z is positive semidefinit
for i = 1:length(O_c)
    F = F + [O_c(i)-C(i,:)*(Z)*C(i,:)'>=0];
end
F = F + [diag(I_c)-C(length(O_c)+1:end,:)*(Z)*C(length(O_c)+1:end,:)'>=0];


%create objective function
ObjecFcn = trace(C*Z*C');
ops=sdpsettings('debug',1,'verbose',1);
sol = optimize(F,ObjecFcn,ops);

%solve for observer gain
% Y = double(Y);
% L = Y*M'*inv(V);

%solve for controller gain
Z = double(Z);
G = double(G);
K = G/Z;

%for i = 1:length(O_c)
%    C(i,:)*Z*C(i,:)'
%end
max(real(eig(A+B*G/Z)))
min(real(eig(A+B*G/Z)))