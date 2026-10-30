function scale = DeformationScale(x,y,ux,uy)
%DEFORMATIONSCALE Scale factor of the displacements in plots of the deformed mesh.

% Scale factor of the displacements for the plot of the deformed mesh
lx = max(x) - min(x);
ly = max(y) - min(y);
scale = min([round(lx/max(abs(ux))) round(ly/max(abs(uy)))])/10;
if isinf(scale) || isnan(scale)
    scale = 1;
end
end
