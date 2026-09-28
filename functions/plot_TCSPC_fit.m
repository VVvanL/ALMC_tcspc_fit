function h = plot_TCSPC_fit(t_bin, TCSPC_trace, r_fitirf, r_fitirf_fit, irf_fit, fit_type)

    h = figure; hold on; grid on

    % plot data trace and fit, IRF fit
    semilogy(t_bin, TCSPC_trace, '.', 'DisplayName','data');
    semilogy(t_bin,r_fitirf_fit,'DisplayName','fit');
    semilogy(t_bin,irf_fit./max(irf_fit)*max(r_fitirf_fit),'DisplayName','irf');
    
    ylim([min([min(r_fitirf_fit) min(TCSPC_trace)]) max(r_fitirf_fit)*1.05]);
    yscale log

    legend;
    xlabel('time (ns)'); ylabel('counts');

    % print fit parameters
    if fit_type == 1

        ant_str = ([char(964), '_1: ',num2str(r_fitirf.taus(1), '%.2f'),' ns ', char(177),' ', num2str(r_fitirf.err_vals.taus(1), '%.2f'),' ns', ...
            newline, 'background: ',num2str(r_fitirf.background, '%.2f'), char(177),' ', num2str(r_fitirf.err_vals.background, '%.2f')]);
        dim = [.66 .7, .1 .1];
        annotation('textbox', 'Position',dim, String = ant_str, FontSize = 13)

    elseif fit_type == 2

        ant_str = ([char(964), '_1: ',num2str(r_fitirf.taus(1), '%.2f'),' ns ', char(177),' ', num2str(r_fitirf.err_vals.taus(1), '%.2f'),' ns', ...
        newline, char(964), '_2: ',num2str(r_fitirf.taus(2), '%.2f'),' ns', char(177),' ', num2str(r_fitirf.err_vals.taus(2), '%.2f'),' ns', ...
        newline, 'background: ',num2str(r_fitirf.background, '%.2f'), char(177),' ', num2str(r_fitirf.err_vals.background,'%.2f')]);
        dim = [.66 .7, .1 .1];
        annotation('textbox', 'Position',dim, String = ant_str, FontSize = 13)   

    end