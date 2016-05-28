odefun = 'Voxel';

L = 0.2;
R_o = 0.005;
th = 0.0025;
rho = 2266;
m_n = 0.11498;

m = 12*L*pi*(R_o^2-(R_o-th)^2)*rho+6*m_n;
k = 218778.999239526;
c = 1;

parameters = {'mass',m;'stiffness',k;'damping',c};

fcn_type = 'c';

init_sys = idgrey(odefun,parameters,fcn_type);

init_sys.Structure.Parameters(1).Free = false;
init_sys.Structure.Parameters(2).Free = true;
init_sys.Structure.Parameters(3).Free = true;

init_sys.Structure.Parameters(1).Minimum = 0;
init_sys.Structure.Parameters(2).Minimum = 0;
init_sys.Structure.Parameters(3).Minimum = 0;

x0 = adjustedTime1Top(1,2:3)';

u = zeros(length(adjustedTime1Top),1);
y = adjustedTime1Top(:,2);
data = iddata(y,u,Ts);

sys = greyest(data,init_sys);

opt = compareOptions('InitialCondition','zero');
compare(data,sys,Inf,opt)