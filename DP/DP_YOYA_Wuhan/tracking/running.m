clc
clear all
v0 = 68.89;
v_max = 90;
v_min = 50;
v_ref = 75;
d = 0.1;
S = 14600;
h = 10;
column=round(S/h)+1;%���ž���ָ�ж�����
alpha = 0.004;% �����źű仯Ȩ��
[X_opt,Z_opt,fval,p,b,mf_sum,a,time_sum,sped_ref]= DP1(v0,v_max,v_min,v_ref,d,S,h,alpha);
save whecow0.7alpha0.004sin.mat
%% plot
figure
width=600;%���bai��������
height=400;%�߶�
left=100;%����Ļ����du��ˮƽ��zhi��
bottem=100;%����Ļ���½Ǵ�ֱ��dao��
set(gcf,'position',[left,bottem,width,height]);
load whecow0.7alpha0.004sin.mat
subplot(3,2,1);
plot(1:10:(column)*10,sped_ref,'-.b','linewidth',2);
hold on
plot(1:10:(column)*10,X_opt,'r','linewidth',2);
xlabel('Distance [m]','Fontname','Times New Roman','fontsize',16);
ylabel('Velocity [m/s]','Fontname','Times New Roman','fontsize',16);
set(gca,'FontSize',16);
xlim([0,14600]);

load whecow0.7alpha0.004sin.mat
subplot(3,2,2);plot(1:10:(column)*10,Z_opt,'r','linewidth',2);
xlabel('Distance [m]','Fontname','Times New Roman','fontsize',16);
ylabel('Fuel Cut Off','Fontname','Times New Roman','fontsize',16); 
set(gca,'FontSize',16);
xlim([0,14600]);
legend('0 - Fuel cut off signal on, 1 - Fuel cut off signal off');

subplot(3,2,3);
load whecow0.7alpha0.004sin.mat
plot(1:10:(column-1)*10,p,'r','linewidth',2);
xlabel('Distance [m]','Fontname','Times New Roman','fontsize',16);
ylabel('Engine Torque [Nm]','Fontname','Times New Roman','fontsize',16);
set(gca,'FontSize',16);

subplot(3,2,4);
load whecow0.7alpha0.004sin.mat
plot(1:10:(column-1)*10,b,'r','linewidth',2);
xlabel('Distance [m]','Fontname','Times New Roman','fontsize',16);
ylabel('Brake Torque [Nm]','Fontname','Times New Roman','fontsize',16);
set(gca,'FontSize',16);
xlim([0,14600]);


subplot(3,2,5);
load whecow0.7alpha0.004sin.mat
plot(1:10:(column)*10,a,'r','linewidth',2);
xlabel('Distance [m]','Fontname','Times New Roman','fontsize',16);
ylabel('Slope [��]','Fontname','Times New Roman','fontsize',16);
set(gca,'FontSize',16);
xlim([0,14600]);