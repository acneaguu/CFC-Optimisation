%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%README:
%%This function is used to create an animated plot of the fitness progression. 
%%The name of the datastruct is used to load the data for the plot
%%The struct should contain a variable F which is used to plot the fitness
%%and a matrix X which is used to plot F as a function of up to 2
%%optimisation variables. With the input 'vars', the position of the
%%columns of the two variables is specified. If vars is empty, the fitness
%%is plotted as function of the number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function animated_plot_fitness(X,F,vars)
%%Refreshrate of the plot
drawrate = 1/144;   

%%Range of ylim w.r.t. final value of F
yres = 25;        

%%If final fitness is nonzero, use custom range
if F(end) > 0
    Flim = F(end)*[(1-yres*0.01),1+yres*0.4];
else
    Flim = [0 1];
end

%%Open a new figure and animated line
if ishandle(1)
    nfig = get(gcf,'Number');
    figure(nfig+1)
else
    figure(1)
end
h = animatedline('Color','red','LineWidth',2);

%%Check for the number of input vars and adapt the axes accordingly
if nargin>2
    if length(vars) > 1
        x1 = X(:,vars(1));
        x2 = X(:,vars(2));
        axis([min(x1),max(x1),min(x2),max(x2),Flim])
        xlabel(sprintf('var %3.1d',vars(1)));
        ylabel(sprintf('var %3.1d',vars(2)));
        zlabel('Fitness')
        view([45 45])
        drawtype = 2;
    else
        x1 = X(:,vars(1));
        axis([0.9*min(x1),1.1*max(x1),Flim])
        xlabel(sprintf('var %3.1d',vars(1)));
        ylabel('Fitness')
        drawtype = 1;
    end
else
    x1 = 1:length(F);
    ylim(Flim)
    xlabel('iteration')
    ylabel('Fitness')
    drawtype = 1;
end 

%%Plot 
a = tic;
for k = 1:length(x1)
    switch drawtype
        case 2
            addpoints(h,x1(k),x2(k),F(k)); 
        case 1 
            addpoints(h,x1(k),F(k));
    end
    b = toc(a);
    if b > drawrate
      drawnow
      a = tic;
    end
end
drawnow
end