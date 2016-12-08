%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = controlledVoxel(t,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% controlledVoxel: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL IN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:    --------------------										
% Calls:	MATLAB 5.2 std fcns
% Inputs:   t - time
%           x - state values
%               The first group of states 1:(length(x))/2 are the actual
%               states. The second group
%               length(x)/2+1:end-m are the "observed" states for the
%               controller. 
% Outputs:  dx - derivates of input state x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: the normal input u is actually part of the states because the
%           state matrix has been augmented by the control u. This causes
%           the u refered to here to be u = v

% Set Global Vairables
global Aaug Baug subM Daug K subK subL subB subBlocks ATildeaug BTildeaug

% Splint the states
x_a = x(1:length(Aaug));
x_c = x(length(Aaug)+1:end);

observ = zeros(size(x_c));
%temp = (Maug*(x_a-x_c));
u = zeros(length(subK),1);
pos = 1;
for i = 1:length(subL)
    [m,n] = size(subB(i).B);
    [s,~] = size(subM(i).M);
    M = [subM(i).M,zeros(s,n);zeros(n,m),eye(n)];
    temp = M*([x_a(2*subBlocks(i)-1:2*subBlocks(i+1));x_a(end-length(subL)+i)]-[x_c(2*subBlocks(i)-1+2*(i-1):2*subBlocks(i+1)+2*(i-1));x_c(end-length(subL)+i)]);
    observ(2*subBlocks(i)-1+2*(i-1):2*subBlocks(i+1)+2*(i-1)) = subL(i).L(1:end-1,:)*temp;
    observ(end-length(subL)+i) = subL(i).L(end,:)*temp;
    
    u(i) = subK(i).K*[x_c(pos:pos+length(subL(i).L)-2);x_c(end-length(subL)+i)];
    pos = pos+length(subL(i).L-1);
end

% Calculate U and the actual state vector
%u = K*x_c;
u = zeros(size(u));
u = K*x_a;

[~,n] = size(Daug);

if t<0.0001
    w = [zeros(n-1,1);1];
else
    w = zeros(n,1);
end

dx_a = Aaug*x_a+Baug*u+Daug*w;

% Calcuate the observed/estimated state vector
dx_c = ATildeaug*x_c + BTildeaug*u + observ;

% Combine derivates
dx = [dx_a;dx_c];