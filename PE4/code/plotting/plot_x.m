function p = plot_x(type, figNr, x, X, params, label_legend, label_title)
%PLOT X Plots for the system states
%   State-state plots state plot and state-time the time plot of the states

    %%% Parse input arguments %%%
    switch nargin
        case 3
            X = [];
            params = {};
            label_legend = '';
            label_title = '';
            
        case 4
            params = {};
            label_legend = '';
            label_title = '';

        case 5
            label_legend = '';
            label_title = '';
        
        case 6
            label_title = '';
            
        case 7
            
        otherwise
            error('Wrong number of inputs!')
    end
    %%%%%%%%%%%%%%%%%%%
    
    nrSteps = size(x,1) - 1;
    switch type
        case 'state-time'
            figure(figNr)
            subplot(2,1,1);
            p = plot(squeeze(x(:,1,:)*180/pi),'color',[params.color,params.alpha],'linewidth',params.lw);
            
            if ~isempty(X)
                hold on;
                plot(linspace(-1, nrSteps+2,2),max(X.V(:,1))*ones(1,2)*180/pi, 'color', 'k','linewidth',2)
                hold on;
                plot(linspace(-1, nrSteps+2,2),min(X.V(:,1))*ones(1,2)*180/pi, 'color', 'k','linewidth',2)
            end
            xlabel('time')
            ylabel('position [deg]')
            xlim([0,nrSteps+2])
            grid()
            if ~isempty(label_legend)
                legend(p, label_legend);
            end
            if ~isempty(label_title)
                title(label_title);
            end
            subplot(2,1,2);
            plot(squeeze(x(:,2,:)*180/pi),'color',[params.color,params.alpha],'linewidth',params.lw)
            if ~isempty(X)
                hold on;
                plot(linspace(-1, nrSteps+2,2),max(X.V(:,2))*ones(1,2)*180/pi, 'color', 'k','linewidth',2)
                hold on;
                plot(linspace(-1, nrSteps+2,2),min(X.V(:,2))*ones(1,2)*180/pi, 'color', 'k','linewidth',2)
            end
            xlabel('time')
            ylabel('velocity [deg/s]')
            xlim([0,nrSteps+2])
            grid()
            set(gcf,'position',[100,100,params.width,2*params.height],'color','white')
        case 'state-state'
            figure(figNr)
            p = plot(squeeze(x(:,1,:)*180/pi), squeeze(x(:,2,:)*180/pi),'color',[params.color,params.alpha],'linewidth',params.lw);
            if ~isempty(X)
                hold on
                X = X*(180/pi);
                X.plot('wire', true, 'linestyle', '-', 'linewidth', 2)
            end
            xlabel('position [deg]')
            ylabel('velocity [deg/s]')
            if ~isempty(label_legend)
                legend(p, label_legend);
            end
            if ~isempty(label_title)
                title(label_title);
            end
            set(gcf,'position',[100,100,params.width,1.5*params.height],'color','white')
        otherwise
            error("invalid plot type!");
    end
end

