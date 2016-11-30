%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = transformationMatrix(ControlLocs,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transformationMatrix: Create transformation matrix for decentralized control	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functions creates a transformation matrix based off of the control
% locations by spliting the states between the control points in half
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:    --------------------										
% Calls:	MATLAB 5.2 std fcns
% Inputs:   ControlLocs - location in the voxel array of control voxels
%           num - the length of the voxel array
% Outputs:  T - transformation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subBlocks = diff([0,diff(ControlLocs)/2+ControlLocs(1:end-1),num]);
for i = 1:(length(subBlocks)-1)
    subBlocks = [subBlocks(1:(2*i-1)),1,subBlocks((2*i):end)];
end

T = [];

for i = 1:length(subBlocks)
    tempT = [];
    if mod(i,2) == 1;
        if i == length(subBlocks)
            for j = 1:length(subBlocks)
                if mod(j,2) == 1
                    if j == length(subBlocks)
                        if i == j
                            tempT = [tempT,eye(2*(subBlocks(i)),2*subBlocks(j))];
                        else
                            tempT = [tempT,zeros(2*(subBlocks(i)),2*subBlocks(j))];
                        end
                    else
                        if i == j
                            tempT = [tempT,eye(2*(subBlocks(i)),2*(subBlocks(j)-1))];
                        else
                            tempT = [tempT,zeros(2*(subBlocks(i)),2*(subBlocks(j)-1))];
                        end
                    end
                else
                    tempT = [tempT,zeros(2*(subBlocks(i)),2)];
                end
            end
        else
            for j = 1:length(subBlocks)
                if mod(j,2) == 1
                    if j == length(subBlocks)
                        if i == j
                            tempT = [tempT,eye(2*(subBlocks(i)-1),2*subBlocks(j))];
                        else
                            tempT = [tempT,zeros(2*(subBlocks(i)-1),2*subBlocks(j))];
                        end
                    else
                        if i == j
                            tempT = [tempT,eye(2*(subBlocks(i)-1),2*(subBlocks(j)-1))];
                        else
                            tempT = [tempT,zeros(2*(subBlocks(i)-1),2*(subBlocks(j)-1))];
                        end
                    end
                else
                    tempT = [tempT,zeros(2*(subBlocks(i)-1),2)];
                end
            end
        end
    else
        tempT2 = [];
        for j = 1:length(subBlocks)
            if mod(j,2) == 1
                if j == length(subBlocks)
                    tempT2 = [tempT2,zeros(2,2*subBlocks(j))];
                else
                    tempT2 = [tempT2,zeros(2,2*(subBlocks(j)-1))];
                end
            else
                if i == j
                    tempT2 = [tempT2,eye(2,2)];
                else
                    tempT2 = [tempT2,zeros(2,2)];
                end
            end
        end
        for j = 1:length(subBlocks)
            if mod(j,2) == 1
                if j == length(subBlocks)
                    tempT = [tempT,zeros(2,2*subBlocks(j))];
                else
                    tempT = [tempT,zeros(2,2*(subBlocks(j)-1))];
                end
            else
                if i == j
                    tempT = [tempT,eye(2,2)];
                else
                    tempT = [tempT,zeros(2,2)];
                end
            end
        end
        tempT = [tempT;tempT2];
    end
    T = [T;tempT];
end