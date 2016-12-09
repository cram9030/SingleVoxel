%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [subK,subG,subZ,K,Ktilde] = DecentralizedLMI(subA,subD,subB,subH,subM,subW,subV,subQ,subR,subI_c,subO_c)
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
Ktilde = zeros(length(subA),tildeLength+length(subA));
%Ltilde = zeros(tildeLength+length(subA),3*length(subA));

%Iterate through each LMI design criteria
pos = 1;
tildePos = 1;
for i = 1:length(subA)
    [subK(i).K,subG(i).G,subZ(i).Z] = LMIDesign(subA(i).A,subD(i).D,subB(i).B,subH(i).H,subM(i).M,subW(i).W,subV(i).V,subQ(i).Q,subR(i).R,subI_c(i).I_c,subO_c(i).O_c);
    K(i,pos:length(subK(i).K)+pos-2) = subK(i).K(1:end-1);
    K(i,end-length(subA)+i) = subK(i).K(end);
    
    Ktilde(i,tildePos:tildePos+length(subK(i).K)-2) = subK(i).K(1:end-1);
    Ktilde(i,end-length(subA)+i) = subK(i).K(end);
    
    pos = pos + length(subK(i).K)-3;
    tildePos = tildePos+length(subK(i).K)-1;
end