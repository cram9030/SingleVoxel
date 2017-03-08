%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = uncontrolledVoxelAero(t,x)
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
global Aaug Baug Daug K

% Calculate U and the actual state vector
%u = K*x_c;
u = K*x;
u = zeros(size(u));

[~,n] = size(Daug);
L = 0:0.283:56.6;
dL = diff(L);
U_inf = 4.02336*rand;
cord = 0.283/.12;
CL = 269.303850944448e-003;
gamma_0 = 2*CL*U_inf*cord/pi;
rho_m = 1.225;
w = zeros(n,1);

for i = 1:200
    w(i) = 1+(1/114*sqrt(3249-(L(i+1)-dL(i)/2)^2)*(L(i+1)-dL(i)/2)+57/2*asin((L(i)-dL(i)/2)/57)-(1/114*sqrt(3249-(L(i)+dL(i)/2)^2)*(L(i)+dL(i)/2)+57/2*asin((L(i)+dL(i)/2)/57)));
end

w = rho_m*U_inf*w;

dx = Aaug*x+Baug*u+Daug*w;