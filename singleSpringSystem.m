%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,B_dist,B,C] = singleSpringSystem(m,k,c,num,SensorLocs,ControlLocs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% singleSpringSystem: Create spring mass damper system of voxels		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function creates the the state, input, and output matrices of an space
% structure made from cuboct vocels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:    --------------------										
% Calls:	MATLAB 5.2 std fcns
% Inputs:   m - mass
%           k - stiffness
%           c- damping
%           ControlLocs - Array of control input locations
%           SensorLocs - Array of sensor locations
%           num - number of voxels
% Outputs:  A - State matrix of size 2*num
%           B_dist - Disturbance input matrix of size 2*num by
%                   length(DisturbLocs)
%           B_cont - Control input matrix of size 2*num by
%                   length(ControlLocs)
%           C - Output matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initalize Matrices
AK = zeros(num);
AC = zeros(num);
B = zeros(2*num,length(ControlLocs));
C = zeros(length(SensorLocs),2*num);

%Populate A matrix
for i = 1:num
    AK(i,i) = -2*k/m;
    AC(i,i) = -2*c/m;
    if i>1
        AK(i,i-1) = k/m;
        AC(i,i-1) = c/m;
    end
    if i<num
        AK(i,i+1) = k/m;
        AC(i,i+1) = c/m;
    end
end
A = [zeros(size(AK)),eye(size(AK));AK,AC];

%Populate disturbance matrix
B_dist = [zeros(num);1/m*eye(num)];

%Populate control input matrices
for i = 1:length(ControlLocs)
    B(num+ControlLocs(i),i) = 1/m;
end

%Create output matrix
for i = 1:length(SensorLocs)
    C(i,SensorLocs(i)) = 1;
end