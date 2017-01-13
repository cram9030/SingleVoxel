N = length(index);
F(N) = struct('cdata',[],'colormap',[]);

figure;
axis([0 199*0.283 -1e-7 0.2e-6])
title('Comparison of Active Voxel Configurations')

ax = gca;
ax.NextPlot = 'replaceChildren';

for j = 1:N
    plot(0:0.283:199*0.283,yUC300(index(j),1:2:400),0:0.283:199*0.283,y25(index(j),1:2:400),0:0.283:199*0.283,yC(index(j),1:2:400),0:0.283:199*0.283,y50(index(j),1:2:400))
    drawnow
    F(j) = getframe(gcf);
end

movie2avi(F,'ActiveVoxelMovie3.avi','compression','None')