function [A,B,C,D] = voxelModel(m,k,c)

A = [-k(2)/c(2),0,0;0,0,1;k(1)/m-c(1)*k(2)/(m*c(2)),-k(1)/m,-c(1)/m];
B = [0;0;1/m];
C = [0,0,1];
D = zeros(1);

end