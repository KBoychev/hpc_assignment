clear all;
close all;
clc;


%% Task 1
% -------------------------------------------------------------

E=210000000000;
I=0.0000144;
L=10;
qy=1000;
Fy=1000;

ts=csvread('task_static.log');
tsa(:,1)=linspace(ts(1,1),ts(end,1),100);
tsa(:,2)=(Fy*tsa(:,1).^2)/(48*E*I).*(3*L-4*tsa(:,1)) + (qy*tsa(:,1).^2)/(24*E*I).*(L-tsa(:,1)).^2;

line_style={'-','--',':','-.','-','--',':','-.','-','--',':','-.'};
line_color=[28,69,135]/255;

figure(1);
hold on;
plot(ts(:,1),ts(:,2),'LineStyle',line_style{1},'Color',line_color,'DisplayName','u_{2} Numerical');
plot(tsa(:,1),tsa(:,2),'LineStyle',line_style{2},'Color',line_color,'DisplayName','u_{2} Analytical');
grid on;
xlabel('L (m)');
ylabel('u_{2} (m)');
legend('show');

%% Task 2
% -------------------------------------------------------------

cndwt=csvread('center_node_defletion_wtr_time.log');
% cndwt_1=csvread('center_node_defletion_wtr_time_1.log');
% cndwt_2=csvread('center_node_defletion_wtr_time_2.log');
% cndwt_3=csvread('center_node_defletion_wtr_time_3.log');

figure(2);
hold on;
plot(cndwt(:,1),cndwt(:,2),'LineStyle',line_style{1},'Color',line_color);
% plot(cndwt_1(:,1),cndwt_1(:,2),'LineStyle',line_style{1},'Color',line_color);
% plot(cndwt_2(:,1),cndwt_2(:,2),'LineStyle',line_style{1},'Color',line_color);
% plot(cndwt_3(:,1),cndwt_3(:,2),'LineStyle',line_style{1},'Color',line_color);
grid on;
ylabel('u_{2} at L/2 (m)');
xlabel('t (s)');


% tde=csvread('task_dynamic_explicit.log');
% tdep=csvread('task_dynamic_explicit_parallel.log');


% tdi=csvread('task_dynamic_explicit.log');
% tdip=csvread('task_dynamic_explicit_parallel.log');