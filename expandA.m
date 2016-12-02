%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Atilde,subA,Btilde,subB,Dtilde,subD,Mtilde,subM,Htilde,subH] = expandA(A,B,D,ControlLocs,num)
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
subB(1).B = [];
subD(1).D = [];
subM(1).M = [];
subH(1).H = [];
Atilde = zeros(length(A)+2*diff(length(ControlLocs)),length(A)+2*diff(length(ControlLocs)));
Btilde = zeros(length(A)+2*diff(length(ControlLocs)),length(ControlLocs));
Dtilde = zeros(length(A)+2*diff(length(ControlLocs)),length(A)/2);
subBlocks = [1,cumsum(diff([0,diff(ControlLocs)/2+ControlLocs(1:end-1),num]))];

%Iterate through all subblocks
for i = 1:length(subBlocks)-1
    %Set primary diagional subBlocks to corrisponded expanded A locations
    Atilde(2*subBlocks(i)-1+2*(i-1):2*subBlocks(i+1)+2*(i-1),2*subBlocks(i)-1+2*(i-1):2*subBlocks(i+1)+2*(i-1)) = A(2*subBlocks(i)-1:2*subBlocks(i+1),2*subBlocks(i)-1:2*subBlocks(i+1));
    Btilde(2*subBlocks(i)-1+2*(i-1):2*subBlocks(i+1)+2*(i-1),:) = B(2*subBlocks(i)-1:2*subBlocks(i+1),:);
    Dtilde(2*subBlocks(i)-1+2*(i-1):2*subBlocks(i+1)+2*(i-1),:) = D(2*subBlocks(i)-1:2*subBlocks(i+1),:);
    
    %Save primary diagional subblocks into structure
    subA(i).A = A(2*subBlocks(i)-1:2*subBlocks(i+1),2*subBlocks(i)-1:2*subBlocks(i+1));
    subB(i).B = B(2*subBlocks(i)-1:2*subBlocks(i+1),i);
    subD(i).D = D(2*subBlocks(i)-1:2*subBlocks(i+1),subBlocks(i):subBlocks(i+1));
    subM(i).M = double(logical(subB(i).B'));
    subH(i).H = double(logical(subB(i).B'));
    
    %Set constraint controls on overlapping subblocks
    if i~=1 && i~=length(subBlocks)-1
        subH(i).H = [subH(i).H;[eye(2),zeros(2,length(subH(i).H)-2)];[zeros(2,length(subH(i).H)-2),eye(2)]];
    elseif i == 1
        subH(i).H = [subH(i).H;[zeros(2,length(subH(i).H)-2),eye(2)]];
    else
        subH(i).H = [subH(i).H;[eye(2),zeros(2,length(subH(i).H)-2)]];
    end
    
    %Insert off diagional subblocks to expanded A
    if i ~= length(subBlocks)-1
        Atilde(2*subBlocks(i+1)+2*(i-1)-1:2*subBlocks(i+1)+2*(i-1),2*subBlocks(i+1)+2*(i-1)+3:2*subBlocks(i+1)+2*(i-1)+4) = A(2*subBlocks(i+1)-1:2*subBlocks(i+1),2*subBlocks(i+1)+1:2*subBlocks(i+1)+2);
        Atilde(2*subBlocks(i+1)+1+2*(i-1):2*subBlocks(i+1)+2*(i-1)+2,2*subBlocks(i+1)+2*(i-1)-3:2*subBlocks(i+1)+2*(i-1)-2) = A(2*subBlocks(i+1)+1:2*subBlocks(i+1)+2,2*subBlocks(i+1)+-1:2*subBlocks(i+1));
    end
end

Mtilde = double(logical(Btilde'));
Htilde = double(logical(Btilde'));