clear
% load('dist.mat');%sampling100
% load('slope.mat');
% load('utopia.mat');
load('sample1.mat');
load('diff1.mat');
dist_samp = 1:1:25000;
cd       = 0.373;%风阻系数
Af       = 2.58;%有效迎风面积
ro       = 1.205;%空气密度
M        = 1870;%整车质量
af       = cd*Af*ro/2/M;%风阻参数
f        = 0.011;%滚动阻力系数
g        = 9.8;%重力加速度
itat     = 0.94;%机械效率
rw       = 0.364;%车轮半径
Ig       = 3.66;
p1       = 0.00385165576001348;
p2       = 1.19220483871151E-6;
p3       = 2.95706278924791E-10;
p4       = 5.45299140436618E-5;
p5       = -1.59515599116976E-9;
p6       = 1.74678055479864E-7; %z = p1*x+p2*x^2+p3*x^3+p4*x*y+p5*x^2*y+p6*x*y^2
%% create grid
clear grd

% Vehicle velocity
grd.Nx{1}    = 101; 
grd.Xn{1}.hi = 25; 
grd.Xn{1}.lo = 15;

% Engine Torque
grd.Nu{1}    = 121;
grd.Un{1}.hi = 120;
grd.Un{1}.lo = 0;

% Gear
grd.Nu{2}    = 1;
grd.Un{2}.hi = 5;
grd.Un{2}.lo = 5;
    
% set initial state
grd.X0{1} = 20;    % initial SOC

% final state constraints
grd.XN{1}.hi = 20.1;
grd.XN{1}.lo = 19.9;

% define problem
clear prb
prb.W{1} = interp1(sample1,diff1,dist_samp*100+100000-1,'previous');%slope
% prb.W{1} = slope_all;
% Utopia = 0.4910*ones(1,18);
% prb.W{2} = interp1(dist800,Utopia800(1,:), dist_samp,'previous'); %Utopia
% prb.W{2} = interp1(dist800,Utopia, dist_samp,'previous'); %Utopia
prb.Ts   = 1; %discrete step size:5m
prb.N    = 25000/prb.Ts ;%(pwr_demand{cycle_num})

% set options
options = dpm();
options.Waitbar = 'on';
options.MyInf = 1e9;
% options.InputType = 'dd';
options.BoundaryMethod = 'none'; % also possible: 'none' or 'LevelSet'
%     if strcmp(options.BoundaryMethod,'Line') 
%        %these options are only needed if 'Line' is used
%        options.Iter = 5;
%        options.Tol = 1e-8;
%        options.FixedGrid = 0;
%     end
tic
[res, dyn] = dpm(@Model,[],grd,prb,options);
toc
save('dynamic_program1v2.mat');
%     cycle_num_str = 'd4h';%num2str(cycle_num);
%     pack_num_str = num2str(par.num_s*par.num_p);
%     filename = datetime('now','Format','MM_dd_yyyy_HHmm');
%     file_name = sprintf('d%s_%s_%s.mat',cycle_num_str,pack_num_str,filename);
%     save(file_name,'dyn2','res2');
%     Fuel = trapz(res2.FFR)/Fuel_Density*Litre_Gallon/3.6;
%     Dist = trapz(speed_vector)*Meter_Mile;
%     res2.MPG = Dist/Fuel;
%     Dist/Fuel
%         
%     str = ['res',DC{k},'=res2;'];
%     eval(str)
%     filename = ['res',DC{k},'_hSC_bs',num2str(par.battery_scale*1000),'_new_sungear_35'];
%     varname = ['res',DC{k}];
%     save([direct filename],varname) 
% figure
% 
% k = 1;
% [min1,minfn1] = fminbnd(@(x)pickindex(x,k),10,20);
% k = 2;
% [min2,minfn2] = fminbnd(@(x)pickindex(x,k),10,20);
% goal = [minfn1,minfn2];
% nf = 2; % number of objective functions
% N = 500; % number of points for plotting
% onen = 1/N;
% x = zeros(N+1,1);
% f = zeros(N+1,nf);
% fun = @simple_mult;
% x0 = 5;
% options = optimoptions('fgoalattain','Display','off');
% for r = 0:N
%     t = onen*r; % 0 through 1
%     weight = [t,1-t];
%     [x(r+1,:),f(r+1,:)] = fgoalattain(fun,x0,goal,weight,...
%         [],[],[],[],[],[],[],options);
% end
% hold on
% %plot(f(:,1),f(:,2),'k.');
% plot(f(:,1),f(:,2),'k.','MarkerFaceColor','g','MarkerEdgeColor','g','LineWidth',2);
% % xlim([0.4,0.8]);
% % ylim([20,70]);
% % xlabel('f_1')
% % ylabel('f_2')
% min_dis = 100000000;
% j2 = 1;
% for i = 1:501
%     dis = sqrt((f(i,1)-f(1,1))^2 + (f(i,2)-f(501,2))^2);
%     if dis<min_dis
%         min_dis = dis;
%         j2 = i;
%     end
% end
% hold on;
% plot(f(1,1),f(501,2),'k.','MarkerSize',20,'MarkerFaceColor','g','MarkerEdgeColor','g','LineWidth',1.5);
% hold on;
% plot(f(j2,1),f(j2,2),'k.','MarkerSize',26,'MarkerFaceColor','g','MarkerEdgeColor','g','LineWidth',1.5);
% hold on
% text(f(j2,1),f(j2,2),'Knee Point','FontSize',24)
% xlabel('Fuel consumption [g/m]','Fontname','Times New Roman','fontsize',35);
% ylabel('Trip time [s]','Fontname','Times New Roman','fontsize',35);
% hold on
% dis      = 5;
% ff(:,1)   = (p1.*kk + p2.*kk.^2.*res.x1(m) + p3.*kk.^3.*res.x1(m).^2 + p4.*kk.*res.u1(m) + p5.*kk.^2.*res.x1(m).*res.u1(m) + p6.*kk.*res.u1(m).^2);
% ff(:,2)   = dis./res.x1(m);
% plot(ff(m,1),ff(m,2),'ro','MarkerSize',15,'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1);
% v = VideoWriter('zero11.avi');
% open(v);
% for m = 1:2721
%     ff(m,1)   = (p1.*kk + p2.*kk.^2.*res.x1(m) + p3.*kk.^3.*res.x1(m).^2 + p4.*kk.*res.u1(m) + p5.*kk.^2.*res.x1(m).*res.u1(m) + p6.*kk.*res.u1(m).^2);
%     ff(m,2)   = dis./res.x1(m);
%     plot(ff(m,1),ff(m,2),'ro','MarkerSize',15,'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1);
%     hold on
%     M2(m) = getframe(gcf);
% end
% writeVideo(v,M2);
% close(v);
% Utopia = zeros(1,1501);
% for i=1:554
%     Utopia(i) = 1.0508;
% end
% for i=555:758
%     Utopia(i) = 1.1163;
% end
% for i=759:888
%     Utopia(i) = 1.0082;
% end
% for i=889:1501
%     Utopia(i) = 1.0508;
% end
% Anti_Utopia = zeros(1,1501);
% for i=1:554
%     Anti_Utopia(i) = 1.5867;
% end
% for i=555:758
%     Anti_Utopia(i) = 1.6572;
% end
% for i=759:888
%     Anti_Utopia(i) = 1.5409;
% end
% for i=889:1501
%     Anti_Utopia(i) = 1.5867;
% end
% figure
% subplot(5,1,1)
% stairs(dist_partial,0.04*ones(1,1501),'Linewidth',2);
% hold on
% stairs(dist_partial,Utopia(1,:),'Linewidth',2);
% set(gca,'FontSize',16);
% legend('UP of Trip Time [s/m]','UP of Fuel Consumption [g/m]')
% % ylabel(hAx(1),'UP of Trip Time [s/m]','Fontname','Times New Roman','fontsize',18)
% % ylabel(hAx(2),'UP of Fuel Consumption [g/m]','Fontname','Times New Roman','fontsize',18);
% subplot(5,1,2)
% stairs(dist_partial,cumtrapz(dist_partial,slope_partial),'Linewidth',2);
% ylabel('Altitude [m]','Fontname','Times New Roman','fontsize',16);
% set(gca,'FontSize',16);
% subplot(5,1,3)
% stairs(dist_partial,res.u1,'Linewidth',2);
% ylabel('Torque [Nm]','Fontname','Times New Roman','fontsize',16);
% set(gca,'FontSize',16);
% subplot(5,1,4)
% stairs(dist_partial,res.u2,'Linewidth',2);
% ylabel('Gear','Fontname','Times New Roman','fontsize',16);
% set(gca,'FontSize',16);
% subplot(5,1,5)
% stairs(dist_partial,res.x1,'Linewidth',2);
% xlabel('Distance [m]','Fontname','Times New Roman','fontsize',16);
% ylabel('Velocity [m/s]','Fontname','Times New Roman','fontsize',16);
% set(gca,'FontSize',16);
% figure
% DP_fuel = [2025.1 2336.4 1899.1 1840.8 1822];
% DP_time = [71.9035 64.2428 77.8785 81.9048 83.4849];
% xx=linspace(68,87);
% yy=spline(DP_time,DP_fuel,xx);
% plot(xx,yy,'r',DP_time,DP_fuel,'o');
% % plot(DP_time,DP_fuel,'Linewidth',2,'LineSmoothing', 'on');
% ylabel('Fuel Consumption [g/m]','Fontname','Times New Roman','fontsize',16);
% xlabel('Trip Time [s/m]','Fontname','Times New Roman','fontsize',16);
% hold on
% plot(76.5218,1914.7 ,'bo'); % Utopia tracking
% hold on
% time = [79.1668, 78.5999, 77.8434, 73.0107];
% fuel = [1845.1, 1858.7, 1877.2, 2011.8];
% plot(time,fuel ,'ro'); % Weighted sum method
% hold on
% plot(80.6565,1865.8 ,'go'); % Utopia tracking
% hold on
% plot(76.7886,1905.2 ,'ro'); % Lyapunov function
% figure;
% xx = 0.1:0.1:8000;
% xx1 = 1:1:8001;
% [hAx,hLine1,hLine2] = plotyy(xx,res.x1,xx1,theta1);
% xlabel('Displacement from origin [m]','Fontname','Times New Roman','fontsize',30)
% ylabel(hAx(1),'Velocity[m/s]','Fontname','Times New Roman','fontsize',30) % left y-axis 
% ylabel(hAx(2),'Slope [ °]','Fontname','Times New Roman','fontsize',30) % right y-axis
% set(gca,'FontSize',38.5);
% % set(hAx(1),'FontSize', 38.5,'Ycolor','k','ylim',[18,21],'ytick',[18:1:21]);
% set(hAx(2),'FontSize', 38.5,'Ycolor','k');
% set(hAx(1),'FontSize', 38.5,'Ycolor','k');
% set(hLine1,'color','r','Linewidth',1);
% set(hLine2,'color','k','Linewidth',1);
% fuel= 0;
% time = 0;
% kk        = 30*Ig/(3.14*rw);
% for m = 1:80000
%     fuelins  = (p1.*kk + p2.*kk.^2.*res.x1(m) + p3.*kk.^3.*res.x1(m).^2 + p4.*kk.*res.u1(m) + p5.*kk.^2.*res.x1(m).*res.u1(m) + p6.*kk.*res.u1(m).^2);
%     timeins  = 1./res.x1(m);
%     fuel = fuel + fuelins*prb.Ts;
%     time = time + timeins*prb.Ts;
% end
% figure
% DP_fuel = [6.615701476275251e+03 6.616909089715974e+03 6.619496617278188e+03 6.700887900004111e+03 7.854300214587453e+03 9.360524956737360e+03];
% DP_time = [5.234408813160238e+02 5.231902612811252e+02 5.225354499640372e+02 5.094579549680124e+02 3.984665203488029e+02 3.268449678840136e+02];
% xx=linspace(3.268449678840136e+02,5.234408813160238e+02);
% yy=spline(DP_time,DP_fuel,xx);
% plot(xx,yy,'r',DP_time,DP_fuel,'o');
% % plot(DP_time,DP_fuel,'Linewidth',2,'LineSmoothing', 'on');
% ylabel('Fuel Consumption [g/m]','Fontname','Times New Roman','fontsize',30);
% xlabel('Trip Time [s/m]','Fontname','Times New Roman','fontsize',30);
% set(gca,'FontSize',38.5);
% hold on
% DP_fuel = [8.032723078479092e+03 7.821223034483986e+03 7.819265253811772e+03];
% DP_time = [3.891273938260799e+02 4.014907381224476e+02 4.016128632630102e+02];
% xx=linspace(3.891273938260799e+02, 4.016128632630102e+02);
% yy=spline(DP_time,DP_fuel,xx);
% plot(xx,yy,'g',DP_time,DP_fuel,'o');