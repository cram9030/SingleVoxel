function [A,B,C,D] = Voxel(m,k,c,Ts)

A = [0 1;-k/m -c/m];
B = [0;1/m];
C = [0,1];
D = zeros(1,1);

end