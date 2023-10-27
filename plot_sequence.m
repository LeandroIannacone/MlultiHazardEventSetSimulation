function [] = plot_sequence(e,resol)

    t_vec = e.t_vec;
    M_vec = e.M_vec;
    typ_vec = e.typ_vec;
    rain_m = e.rain_m;

    figure(1)
    % Mainshocks
    t_ms_vec = t_vec(typ_vec == 1);
    m_ms_vec = M_vec(typ_vec == 1);
    % Aftershocks
    t_as_vec = t_vec(typ_vec == 2);
    m_as_vec = M_vec(typ_vec == 2);
    % Rain
    t_r_vec = t_vec(typ_vec == 4);
    m_r_vec = M_vec(typ_vec == 4);
    % Landslide
    t_l_vec = t_vec(typ_vec == 5);
    
    scatter(t_as_vec/resol,ones(1,length(t_as_vec)),exp(m_as_vec)/2,'magenta')
    hold on
    scatter(t_ms_vec/resol,ones(1,length(t_ms_vec)),exp(m_ms_vec),'b','LineWidth',1.5)
    scatter(t_r_vec/resol,ones(1,length(t_r_vec)),(rain_m(:,1)+1),'k')
    scatter(t_ms_vec/resol,2*ones(1,length(t_ms_vec)),exp(m_ms_vec),'b','LineWidth',1.5)
    scatter(t_as_vec/resol,3*ones(1,length(t_as_vec)),exp(m_as_vec)/2,'magenta')
    scatter(t_r_vec/resol,ones(1,length(t_r_vec)),(rain_m(:,2)+1),'k')
    scatter(t_r_vec/resol,4*ones(1,length(t_r_vec)),(rain_m(:,1)+1),'k')
    scatter(t_r_vec/resol,5*ones(1,length(t_r_vec)),(rain_m(:,2)+1),'k')
    scatter(t_l_vec/resol,ones(1,length(t_l_vec)),'+','r','Linewidth',1.5)
    scatter(t_l_vec/resol,6*ones(1,length(t_l_vec)),'+','r','Linewidth',1.5)
    hold off
    xlim([0 50])
    yticks([1 2 3 4 5 6])
    xlabel('Time [yr]','Interpreter','latex','FontSize',16)
    yticklabels({'All','Main Shock','Aftershock','Rain Intensity','Rain Duration','Landslide'})
    ax = gca;
    ax.FontSize = 16;
    ax.TickLabelInterpreter = "latex";

end