function plot_PE7_3D(figNr, bo, X, params)
%PLOT PE 7 - 3D
    %%% Parse input arguments %%%
    switch nargin
        case 4
                        
        otherwise
            error('Wrong number of inputs!')
    end
    %%%%%%%%%%%%%%%%%%%
    
    figure(figNr); hold on;
    title('Gaussian Process Approximation')
    [y, ~] = bo.sample(X);
    data = bo.get_data();
    [theta, mean, ~] = bo.get_estimate(X);
    
    [Xq, Yq] = meshgrid(params.ctrl.range(1,1):0.1:params.ctrl.range(1,2), params.ctrl.range(2,1):0.1:params.ctrl.range(2,2));
    Eq = griddata(X(1,:), X(2,:), y', Xq, Yq);
    p1 = mesh(Xq, Yq, Eq);
    p2 = plot3(data.x(1,:), data.x(2,:), data.y','b+');   
    p3 = plot3(theta(1,:),theta(2,:), mean, 'rx','MarkerSize',10);
    
    xlabel('P gain'); ylabel('I gain'); zlabel('mean GP value');
    legend([p2,p1,p3],{'Data','Fitted GP Mean','Best Parameter'},'Location','nw');
    set(gcf,'position',[100,100,0.7*params.plot.width,1.4*params.plot.height],'color','white')
end

