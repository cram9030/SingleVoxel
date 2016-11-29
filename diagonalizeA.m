%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ahat,Dhat,Bhat,Hhat,Mhat] = diagonalizeA(A,D,B,H,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diagonalizeA: diagonalize the A matrix procuced by singleSpringSystem		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function diagonalizes A to create decentralized form
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

%Initalize output matrices
Ahat = zeros(size(A));
Dhat = zeros(size(D));
Bhat = zeros(size(B));
Hhat = zeros(size(H));
Mhat = zeros(size(M));

for i = 1:length(A)/2
     Ahat(2*(i-1)+1,2*(i-1)+2) = A(i,length(A)/2+i);
     index = find(A(length(A)/2+i,:));
     temp = A(length(A)/2+i,index);
     if length(index) == 4
         if i == 1
             Ahat(2,1:4) = [temp(1),temp(3),temp(2),temp(4)];
         else
             Ahat(end,length(A)-3:end) = [temp(1),temp(3),temp(2),temp(4)];
         end
     else
         Ahat(2*i,2*(i-2)+1:2*(i-2)+6) = [temp(1),temp(4),temp(2),temp(5),temp(3),temp(6)];
     end
     Dhat(2*i,:) = D(length(A)/2+i,:);
     Bhat(2*i,:) = B(length(A)/2+i,:);
     Hhat(:,2*i) = H(:,length(A)/2+i);
     Hhat(:,2*(i-1)+1) = H(:,i);
     Mhat(:,2*i) = H(:,length(A)/2+i);
     Mhat(:,2*(i-1)+1) = H(:,i);
end