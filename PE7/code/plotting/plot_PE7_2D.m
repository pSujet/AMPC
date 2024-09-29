function plot_PE7_2D(figNr, bo, f, y_hat, x_hat, params)
%PLOT PE 7 - 2D
    %%% Parse input arguments %%%
    switch nargin
        case 6
                        
        otherwise
            error('Wrong number of inputs!')
    end
    %%%%%%%%%%%%%%%%%%%
    
    figure(figNr); hold on;
    title('Gaussian Process Approximation')
    x_plot=0:0.01:1; [mean_plot, std_plot] = bo.sample(x_plot);
    data = bo.get_data();
    [~, idx] = min(y_hat);
    theta = x_hat(idx);
    plot(x_plot, mean_plot+1.96*std_plot, 'k--', 'LineWidth',1.5);
    plot(x_plot, mean_plot-1.96*std_plot, 'k--', 'LineWidth',1.5);
    p1 = fill([x_plot'; flip(x_plot')], [mean_plot+1.96*std_plot; flip(mean_plot-1.96*std_plot)], [220,220,220]/256, 'FaceAlpha', 0.8, 'EdgeAlpha',0);
    p2 = plot(x_plot,mean_plot,'k','LineWidth',1.5);
    p3 = plot(data.x,data.y,'bx','MarkerSize',15);
    p4 = plot(0:0.01:1, f(0:0.01:1),'r');
    if ~isempty(y_hat) && ~isempty(x_hat)
        p5 = scatter(x_hat, y_hat, 8, 'mo');
        p6 = plot(theta,f(theta),'rx','MarkerSize',20);
    end
    xlabel('x'); ylabel('y'); grid();
    if ~isempty(y_hat) && ~isempty(x_hat)
        legend([p3,p2,p1(1),p4,p6,p5],{'Data','Mean','95% Conf. Interval','True Function','Next Point','Acqu. Fn. Samples'},'Location','nw');
    else
        legend([p3,p2,p1(1),p4],{'Data','Mean','95% Conf. Interval','True Function'},'Location','nw');
    end
    set(gcf,'position',[100,100,0.7*params.width,1.4*params.height],'color','white')
end

