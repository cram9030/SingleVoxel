function dx = voxelVib(t,x)

global A

dx = zeros(size(x));
dx = A*x;