function plot_rec08(figNr, data, blr, x, params)
%PLOT plot function for recitation 8

    figure(figNr)
    hold on
    colormap winter;
    
    x_plot=(0:0.1:10);
    
    [mean_plot, var_plot] = blr.predict(x_plot);
    
    plot(x_plot, mean_plot+3*sqrt(var_plot), 'k--', 'LineWidth',1.5);
    plot(x_plot, mean_plot-3*sqrt(var_plot), 'k--', 'LineWidth',1.5);

    p4 = fill([x_plot'; flip(x_plot')], [mean_plot+3*sqrt(var_plot); flip(mean_plot-3*sqrt(var_plot))], [220,220,220]/256, 'FaceAlpha', 0.8, 'EdgeAlpha',0);
    p1 = plot(data.x,data.y,'b.','MarkerSize',25);
    p2 = plot(x_plot, mean_plot,'k','LineWidth',1.5);

    if ~isempty(x)
        p3 = plot(x(1),x(2), 'r.', 'MarkerSize', 30);
        legend([p1,p2,p4,p3],{'Data','Mean','Uncertainty','Prediction'},'Location','nw');
    else
        legend([p1,p2,p4],{'Data','Mean','Uncertainty'},'Location','nw')
    end
    hold off
    xlabel("Chocolate consumption per capita per year");
    ylabel("Nobel laureates per 10 Mio.");
    set(gcf,'position',[0,0, params.width,1.5*params.height],'color','white')
end

