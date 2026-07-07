function plot_trajectory_heights(model)

figure
clf
for ii=1:numel(model.obs)
    plot(model.obs(ii).t_sec,-model.obs(ii).r_ned(:,3),'k.')
    hold on;
end
    plot(model.val_obs.t_sec,-model.val_obs.r_ned(:,3),'r.')
grid minor;
ylabel('Height [m]')
xlabel('Time [s]')


end


