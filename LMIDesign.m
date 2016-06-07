function [L,G,Z] = LMIDesign(A,D,B,H,M,W,V,I_c,O_c)

yalmip('clear');
[m,n] = size(B);
[~,r] = size(D);
A = [A,B;zeros(n,m),zeros(n)];
B = [zeros(size(B));eye(n)];
D = [D;zeros(n,r)];
M = [M,zeros(length(O_c),n)];
C = [H,zeros(length(O_c),length(I_c));zeros(length(I_c),m),eye(length(I_c))];

Y = sdpvar(length(A),length(A));
Z = sdpvar(length(A),length(A));
G = sdpvar(n,m+n);

%estimator
F = [[A*Y+Y*A'+D*W*D' Y*M';M*Y V]<0];
F = F + [Y>0];%Enforce positive semidefinite Y

%Controller
F = F + [[Z*A'+A*Z+G'*B'+B*G  Y*M'*inv(V)*sqrt(V);sqrt(V)'*(Y*M'*inv(V))' -eye(size(V))]<0];
F = F + [Z > 0]; % %Z is positive semidefinit
F = F + [diag(O_c)-C(1:length(O_c),:)*(Y+Z)*C(1:length(O_c),:)'>=0];
F = F + [diag(I_c)-C(length(O_c)+1:end,:)*(Y+Z)*C(length(O_c)+1:end,:)'>=0];

% %estimator
% O = [[A*Y+Y*A'+D*W*D' Y*M';M*Y V]<0];
% O = O + [Y>0];%Enforce positive semidefinite Y
% ObjecFcn = trace(C*Y*C');
% sol = solvesdp(O,ObjecFcn);
% 
% %solve for observer gain
% Y = double(Y);
% L = Y*M'*inv(V);
% 
% %Controller
% F = [[Z*A'+A*Z+G'*B'+B*G  L*sqrt(V);sqrt(V)'*L' -eye(size(V))]<0];
% F = F + [Z > 0]; % %Z is positive semidefinit
% F = F + [diag(O_c)-C(1:length(O_c),:)*(Y+Z)*C(1:length(O_c),:)'>=0];
% F = F + [diag(I_c)-C(length(O_c)+1:end,:)*(Y+Z)*C(length(O_c)+1:end,:)'>=0];

%create objective function
ObjecFcn = trace(C*(Y+Z)*C');
ops=sdpsettings('debug',1,'verbose',1);
sol = solvesdp(F,ObjecFcn,ops);

%solve for observer gain
Y = double(Y);
L = Y*M'*inv(V);

%solve for controller gain
Z = double(Z);
G = double(G);
%K = G*pinv(Z);