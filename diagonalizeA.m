%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Ahat,Dhat,Bhat,Hhat,Mhat] = diagonalizeA(A,D,B,H,M)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diagonalizeA: diagonalize the A matrix procuced by singleSpringSystem		
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function diagonalizes A and shift other matrices to corrispond
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:    --------------------										
% Calls:	MATLAB 5.2 std fcns
% Inputs:   A - State Matrix
%           D - Disturbance Matrix
%           B - Input Matrix
%           H - Control State Matrix
%           M - Output Matrix
% Outputs:  Ahat - Diagonalized State Matrix
%           Dhat - Adjusted Disturbance Matrix
%           Bhat - Adjusted Input Matrix
%           Hhat - Adjusted Control State Matrix
%           Mhat - Adjusted Output Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initalize output matrices
Ahat = zeros(size(A));
Dhat = zeros(size(D));
Bhat = zeros(size(B));
Hhat = zeros(size(H));
Mhat = zeros(size(M));

%Iterate through subblock matrices
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