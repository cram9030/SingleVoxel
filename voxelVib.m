function dx = voxelVib(t,x)

global c k m

dx = zeros(2,1);
dx(1) = x(2);
dx(2) = 1/m*(-k*x(1));%-c*x(2));