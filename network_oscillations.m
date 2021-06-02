% x,y cycling (set by previous model to create oscillations)
a = 0.01;
b = 0.1;
c = 0.01;
d = 0.1;
e = 100;

%transport terms
beta = 0.5; %z diffuses in an unbiased way throughout tree
delta = 0.05;
%delta = [0 0.1 10 1000];

numsets = length(delta);
num_cells = 4; %for other values, the plots will need some tweaking, default is for 4 cell networks
network_type = 'tree'; %valid inputs are 'tree' (must be power of 2 less than 32), 'star', and 'line'
adj_beta = build_adj(beta,num_cells, network_type);

%can build adjacency matrix and create graph for visualization
%adj_mat = 2*adj_beta + diag(sum(abs(2*adj_beta))/2);
%box on; axis square; plot(graph(adj_mat),'layout','force','MarkerSize',16,'LineWidth',2,'EdgeColor','k','NodeFontSize',12);

%set up initial conditions
x_offset = 0.001; %left endpoint of x interval
x_width = 6; %length of x interval
y_offset = 0.001; %lower endpoint for y interval
y_width = 6; %length of y interval
rand_genx = x_width*rand(1,num_cells-1)+x_offset;
rand_geny = y_width*rand(1,num_cells-1)+y_offset;

%used to keep delta scans consistent across trials
%rand_genx = [3.6129    1.5788    3.9255];
%rand_geny = [4.1363    4.4899    2.7042];

%load unperturbed limit cycle points for comparisons
load('unp.mat', 'unpert');
load('unpy.mat', 'unperty');

figure;
for qq = 1:length(delta) %for length of delta, make plots for each value
    
    v = zeros(2*num_cells,1); %set up ICs for time evolution
    v(3:2:end) = rand_genx;
    v(4:2:end) = rand_geny;
    
    %run initial transient and then another round (should be converged if
    %limit cycle by this point)
    [t1,output1] = ode45(@(t,v) simple_cyc_ode(v,adj_beta,a,b,c,d,e,delta(qq)), [0 30000], v);
    [t,output] = ode45(@(t,v) simple_cyc_ode(v,adj_beta,a,b,c,d,e,delta(qq)),[0 10000],output1(end,:));
    
    %convert output for x,y
    x = output(:,1:2:size(output,2));
    y = output(:,2:2:size(output,2));
    
    %plot x and y for each nurse cell in different ways
    %phase plane for x_i,y_i
    subplot(3,numsets,qq);
    set(gca,'fontsize',18)
    box on; grid on; hold on;
    axis square;
    xlabel('$X(t)$','interpreter','latex')
    ylabel('$Y(t)$','interpreter','latex')
    ylim([0 6])
    h = gca;
    h.FontSize = 24;
    title(strcat(['$\delta = $' ' ' num2str(delta(qq))]),'interpreter','latex','FontSize',32)
    plot(x(:,2),y(:,2),'Color',[92/255 157/255 178/255],'LineWidth',3);
    plot(x(:,3),y(:,3),'Color',[180/255 67/255 59/255],'LineWidth',3);
    plot(x(:,4),y(:,4),'Color',[251/255 192/255 52/255],'LineWidth',3);
    plot(unpert,unperty,'k-');
    hold off;
    
    %plot of time evolution of x_i
    subplot(3,numsets,qq+numsets);
    set(gca,'fontsize',18)
    box on; grid on; hold on;
    axis square;
    xlabel('Time, $t$','interpreter','latex')
    ylabel('$X(t)$','interpreter','latex')
    xlim([0 1500])
    ylim([0 6])
    h = gca;
    h.FontSize = 24;
    plot(t,x(:,2),'Color',[92/255 157/255 178/255],'LineWidth',3);
    plot(t,x(:,3),'Color',[180/255 67/255 59/255],'LineWidth',3);
    plot(t,x(:,4),'Color',[251/255 192/255 52/255],'LineWidth',3);
    hold off;
    
    %3D phase plane for all x_i
    subplot(3,numsets,qq+2*numsets);
    plot3(x(:,2),x(:,3),x(:,4),'k-','LineWidth',3);
    set(gca,'fontsize',18)
    hold on; grid on; box on; axis square;
    xlabel('$X_2(t)$','interpreter','latex')
    ylabel('$X_3(t)$','interpreter','latex')
    zlabel('$X_4(t)$','interpreter','latex')
    xticks([0 2 4 6])
    yticks([0 2 4 6])
    zticks([0 2 4 6])
    axis([0 6 0 6 0 6])
    h = gca;
    h.FontSize = 24;
    
end
