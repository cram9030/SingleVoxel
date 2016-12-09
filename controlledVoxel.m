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
global Aaug Baug Daug K

% Calculate U and the actual state vector
%u = K*x_c;
%u = zeros(size(u));
u = K*x;

[~,n] = size(Daug);

if t<0.0001
    w = [zeros(n-1,1);1];
else
    w = zeros(n,1);
end

dx = Aaug*x+Baug*u+Daug*w;