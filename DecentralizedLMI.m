%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [subL,subK,subG,subZ,K,L] = DecentralizedLMI(subA,subD,subB,subH,subM,subW,subV,subI_c,subO_c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DecentralizedLMI: Computes the decralized LMI matrices	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculated the observer and controller gain matrices for
% the decentralized system that is given
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Author:	Nick Cramer, UCSC, Department of Computer Engineering			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called by:    --------------------										
% Calls:	MATLAB 5.2 std fcns
% Inputs:   subA - State Matrix Structure
%           subD - Disturbance Matrix Structure
%           subB - Input Matrix Structure
%           subH - Control State Matrix Structure
%           subM - Output Matrix Structure
%           subW - Disturbance Noise Matrix Structure
%           subV - Sensor Noise Matrix Structure
%           subI_c - Input Covariance Weighting Matrix Structure
%           subO_c - Output Covariance Weighting Matrix Structure
% Outputs:  subL - Observer Gain Structure
%           subK - Controller Gain Structure
%           subG - G Strcture
%           subZ - Z Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Intialize Strcutres
subL(1).L = [];
subG(1).G = [];
subZ(1).Z = [];
subK(1).K = [];

tildeLength = 0;
for i = 1:length(subA)
    tildeLength = tildeLength+length(subA(i).A);
end
K = zeros(length(subA),tildeLength-length(subA)+2);
L = zeros(tildeLength-length(subA)+2,2*length(subA));

%Iterate through each LMI design criteria
pos = 1;
for i = 1:length(subA)
    [subL(i).L,subG(i).G,subZ(i).Z] = LMIDesign(subA(i).A,subD(i).D,subB(i).B,subH(i).H,subM(i).M,subW(i).W,subV(i).V,subI_c(i).I_c,subO_c(i).O_c);
    subK(i).K = subG(i).G/subZ(i).Z;
    K(i,pos:length(subK(i).K)+pos-2) = subK(i).K(1:end-1);
    K(i,end-length(subA)+i) = subK(i).K(end);
    L(pos:length(subL(i).L)+pos-2,i) = subL(i).L(1:end-1,1);
    L(end-length(subA)+i,i) = subL(i).L(end,1);
    L(pos:length(subL(i).L)+pos-2,i+length(subA)) = subL(i).L(1:end-1,2);
    L(end-length(subA)+i,i+length(subA)) = subL(i).L(end,2);
    pos = pos + length(subK(i).K)-3;
end