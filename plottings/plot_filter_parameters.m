function plot_filter_parameters(obs,filter_par)
%PLOT_FILTER_PARAMETERS Plot stored Kalman filter parameters.
%
% Usage:
%   plot_filter_parameters(obs,filter_par)
%
% The function plots
%   1) empirical NIS CDF together with chi-square CDF with 2 DOF,
%   2) temporal bias g with +/- 2 standard deviations,
%   3) sensor bias with +/- 2 standard deviations,
%   4) front orientation bias components with +/- 2 standard deviations,
%   5) back orientation bias components with +/- 2 standard deviations.



    nFlight = numel(obs);
    nOri = 3;

    % NIS empirical CDF and theoretical chi-square CDF with 2 DOF.
    figure()
    clf
for ii = 1:nFlight
    subplot(nFlight,1,ii)

    nis = filter_par(ii).NIS;

    % Empirical CDF
    [F,z] = ecdf(nis);

    plot(z,F,'LineWidth',1.5)
    hold on

    % Theoretical chi-square CDF (2 DOF)
    z_th = linspace(0,max(10,max(z)),500);
    plot(z_th,chi2cdf(z_th,2),'--','LineWidth',1.5)

    hold off
    xlim([0 10])
    ylim([0 1])
    grid minor
    xlabel('NIS')
    ylabel('CDF')
    legend('Estimated','\chi^2(2)','Location','southeast')

    title(sprintf('Flight %d: Mean NIS = %.2f',ii,mean(nis)))
end

    % Temporal bias g.
    figure();    
    clf
    for ii = 1:nFlight
        subplot(nFlight,1,ii)
        t = obs(ii).t_sec(:);
        mu = filter_par(ii).g(:);
        sig = sqrt(max(filter_par(ii).g_var(:),0));

        plot_state_with_2sigma(t,mu,sig)
        grid minor
        xlabel('Time [s]')
        ylabel('g')
        title(sprintf('Flight %d: Temporal bias',ii))
        xlim([t(2) t(end)])
    end

    % Front orientation bias.
    figure()
    clf
    for ii = 1:nFlight
        for jj = 1:nOri
            subplot(nFlight,nOri,(ii-1)*nOri + jj)
            t = obs(ii).t_sec(:);
            mu = filter_par(ii).ori_bias_front(:,jj);
            sig = sqrt(max(filter_par(ii).ori_bias_front_cov_diag(:,jj),0));

            plot_state_with_2sigma(t,mu,sig)
            grid minor
            xlabel('Time [s]')
            ylabel(sprintf('ori front %d',jj))
            title(sprintf('Flight %d: Front ori bias %d',ii,jj))
            ylim([-3e-4 3e-4])
        end
    end

    % Back orientation bias.
    figure()
    clf
    for ii = 1:nFlight
        for jj = 1:nOri
            subplot(nFlight,nOri,(ii-1)*nOri + jj)
            t = obs(ii).t_sec(:);
            mu = filter_par(ii).ori_bias_back(:,jj);
            sig = sqrt(max(filter_par(ii).ori_bias_back_cov_diag(:,jj),0));

            plot_state_with_2sigma(t,mu,sig)
            grid minor
            xlabel('Time [s]')
            ylabel(sprintf('ori back %d',jj))
            title(sprintf('Flight %d: Back ori bias %d',ii,jj))
            ylim([-3e-4 3e-4])
        end
    end
end

function plot_state_with_2sigma(t,mu,sig)
%PLOT_STATE_WITH_2SIGMA Plot state estimate and +/- 2 standard deviations.

    n = min([numel(t),numel(mu),numel(sig)]);
    t = t(1:n);
    mu = mu(1:n);
    sig = sig(1:n);

    plot(t,mu,'r','LineWidth',1.5)
    hold on
    plot(t,mu + 2*sig,'k','LineWidth',1.0)
    plot(t,mu - 2*sig,'k','LineWidth',1.0)
    hold off
end
