%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Atilde,subA] = expandA(A,ControlLocs,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% expandA: Expands A for decentralized control of overlaping information structure	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function expands A for decentralizd control
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:    --------------------										
% Calls:	MATLAB 5.2 std fcns
% Inputs:   A - Banded State Matrix
%           ControlLocs - Control voxel locations
%           num - Total number of voxels
% Outputs:  Atilde - Expanded A matrix
%           subA - array of matrices
%                   A - subBlock matrix to be controlled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Check that Matrix is appropriately banded for a spring mass representation
[lower,upper] = bandwidth(A);
if lower ~= 3 || upper ~= 2
    error('Matrix must be lower bounded 3 and upper bounded 2');
end

%Initialize Parameters
subA(1).A = [];
Atilde = zeros(length(A)+2*diff(length(ControlLocs)),length(A)+2*diff(length(ControlLocs)));
subBlocks = [1,cumsum(diff([0,diff(ControlLocs)/2+ControlLocs(1:end-1),num]))];

%Iterate through all subblocks
for i = 1:length(subBlocks)-1
    %Set primary diagional subBlocks to corrisponded expanded A locations
    Atilde(2*subBlocks(i)-1+2*(i-1):2*subBlocks(i+1)+2*(i-1),2*subBlocks(i)-1+2*(i-1):2*subBlocks(i+1)+2*(i-1)) = A(2*subBlocks(i)-1:2*subBlocks(i+1),2*subBlocks(i)-1:2*subBlocks(i+1));
    %Save primary diagional subblocks into structure
    subA(i).A = A(2*subBlocks(i)-1:2*subBlocks(i+1),2*subBlocks(i)-1:2*subBlocks(i+1));
    
    %Insert off diagional subblocks to expanded A
    if i ~= length(subBlocks)-1
        Atilde(2*subBlocks(i+1)+2*(i-1)-1:2*subBlocks(i+1)+2*(i-1),2*subBlocks(i+1)+2*(i-1)+3:2*subBlocks(i+1)+2*(i-1)+4) = A(2*subBlocks(i+1)-1:2*subBlocks(i+1),2*subBlocks(i+1)+1:2*subBlocks(i+1)+2);
        Atilde(2*subBlocks(i+1)+1+2*(i-1):2*subBlocks(i+1)+2*(i-1)+2,2*subBlocks(i+1)+2*(i-1)-3:2*subBlocks(i+1)+2*(i-1)-2) = A(2*subBlocks(i+1)+1:2*subBlocks(i+1)+2,2*subBlocks(i+1)+-1:2*subBlocks(i+1));
    end
end