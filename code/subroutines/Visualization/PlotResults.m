function RES = PlotResults(model,settings,u)
%PLOTRESULTS Contour plots of a displacement field on the deformed mesh.
%
%   RES = PlotResults(model, settings, u) plots the fields selected in the
%   settings (plot_disp: magnitude, ux, uy; plot_stress: von Mises,
%   sigma_x, sigma_y, tau_xy; plot_strain: epsilon_x, epsilon_y, gamma_xy)
%   with Contour_plot. The displacements are scaled so that the largest
%   one is about a tenth of the model size. Returns the output of
%   Contour_plot, or [] if nothing is plotted.

fields = {};
if settings.plot_disp
    fields = [fields, {'mag','ux','uy'}];
end
if settings.plot_stress
    fields = [fields, {'Mises','sigma_x','sigma_y','tau_xy'}];
end
if settings.plot_strain
    fields = [fields, {'epsilon_x','epsilon_y','gamma_xy'}];
end
if isempty(fields)
    RES = [];
    return
end

scale = DeformationScale(model.xCoor,model.yCoor,u(1:2:end),u(2:2:end));

% Face transparency and nodal averaging of stresses (1: yes, 0: no)
n_alph   = 1;
flag_ave = 0;
RES = Contour_plot(model.p,model.xCoor,model.yCoor,u,scale,fields,...
      model.Topology_Node,model.C,model.material_elset,flag_ave,n_alph);

end
