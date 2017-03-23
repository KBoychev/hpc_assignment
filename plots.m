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

%% Task 1
% -------------------------------------------------------------

figure('position', [100, 100, 700, 225]);
hold on;
plot(ts(:,1),ts(:,2),'LineStyle',line_style{1},'Color',line_color,'DisplayName','Numerical');
plot(tsa(:,1),tsa(:,2),'LineStyle',line_style{2},'Color',line_color,'DisplayName','Analytical');
grid on;
box off;
xlabel('L (m)');
ylabel('u_{2} (m)');
legend('show');

%% Task 2
% -------------------------------------------------------------

e_cndwt=csvread('explicit_center_node_defletion_wtr_time_TleqT.log');
e_cndwt_1=csvread('explicit_center_node_defletion_wtr_time_Tleq05T.log');
e_cndwt_2=csvread('explicit_center_node_defletion_wtr_time_Tleq06T.log');
e_cndwt_3=csvread('explicit_center_node_defletion_wtr_time_Tleq07T.log');
e_cndwt_4=csvread('explicit_center_node_defletion_wtr_time_Tleq08T.log');
e_cndwt_5=csvread('explicit_center_node_defletion_wtr_time_Tleq09T.log');

figure('position', [100, 100, 700, 225]);
hold on;
plot(e_cndwt(:,1),e_cndwt(:,2),'LineStyle',line_style{1},'Color',line_color,'DisplayName','T_{l}=T');
plot(e_cndwt_1(:,1),e_cndwt_1(:,2),'LineStyle',line_style{2},'Color',line_color,'DisplayName','T_{l}=0.5T');
plot(e_cndwt_2(:,1),e_cndwt_2(:,2),'LineStyle',line_style{3},'Color',line_color,'DisplayName','T_{l}=0.6T');
plot(e_cndwt_3(:,1),e_cndwt_3(:,2),'LineStyle',line_style{4},'Color',line_color,'DisplayName','T_{l}=0.7T');
plot(e_cndwt_4(:,1),e_cndwt_4(:,2),'LineStyle',line_style{5},'Color',line_color,'DisplayName','T_{l}=0.8T');
plot(e_cndwt_5(:,1),e_cndwt_5(:,2),'LineStyle',line_style{6},'Color',line_color,'DisplayName','T_{l}=0.9T');
grid on;
box off;
ylabel('u_{2} at L/2 (m)');
xlabel('t (s)');
legend('show');

T=1;
Nt=8000;
dt=T/(Nt-1);

i=round(0.5/dt);
u_amplitude(1)=max(e_cndwt_1(i:end,2))-min(e_cndwt_1(i:end,2));

i=round(0.6/dt);
u_amplitude(2)=max(e_cndwt_2(i:end,2))-min(e_cndwt_2(i:end,2));

i=round(0.7/dt);
u_amplitude(3)=max(e_cndwt_3(i:end,2))-min(e_cndwt_3(i:end,2));

i=round(0.8/dt);
u_amplitude(4)=max(e_cndwt_4(i:end,2))-min(e_cndwt_4(i:end,2));

i=round(0.9/dt);
u_amplitude(5)=max(e_cndwt_5(i:end,2))-min(e_cndwt_5(i:end,2));

Tl=[0.5,0.6,0.7,0.8,0.9];

figure('position', [100, 100, 700, 225]);
loglog(Tl,u_amplitude,'LineStyle',line_style{1},'Color',line_color);
grid on;
box off;

%% Task 3
% -------------------------------------------------------------

i_cndwt=csvread('implicit_center_node_defletion_wtr_time_TleqT.log');
i_cndwt_1=csvread('implicit_center_node_defletion_wtr_time_Tleq05T.log');
i_cndwt_2=csvread('implicit_center_node_defletion_wtr_time_Tleq06T.log');
i_cndwt_3=csvread('implicit_center_node_defletion_wtr_time_Tleq07T.log');
i_cndwt_4=csvread('implicit_center_node_defletion_wtr_time_Tleq08T.log');
i_cndwt_5=csvread('implicit_center_node_defletion_wtr_time_Tleq09T.log');

figure('position', [100, 100, 700, 225]);
hold on;
plot(i_cndwt(:,1),i_cndwt(:,2),'LineStyle',line_style{1},'Color',line_color,'DisplayName','T_{l}=T');
plot(i_cndwt_1(:,1),i_cndwt_1(:,2),'LineStyle',line_style{2},'Color',line_color,'DisplayName','T_{l}=0.5T');
plot(i_cndwt_2(:,1),i_cndwt_2(:,2),'LineStyle',line_style{3},'Color',line_color,'DisplayName','T_{l}=0.6T');
plot(i_cndwt_3(:,1),i_cndwt_3(:,2),'LineStyle',line_style{4},'Color',line_color,'DisplayName','T_{l}=0.7T');
plot(i_cndwt_4(:,1),i_cndwt_4(:,2),'LineStyle',line_style{5},'Color',line_color,'DisplayName','T_{l}=0.8T');
plot(i_cndwt_5(:,1),i_cndwt_5(:,2),'LineStyle',line_style{6},'Color',line_color,'DisplayName','T_{l}=0.9T');
grid on;
box off;
ylabel('u_{2} at L/2 (m)');
xlabel('t (s)');
legend('show');

T=1;
Nt=1000;
dt=T/(Nt-1);

i=round(0.5/dt);
u_amplitude(1)=max(i_cndwt_1(i:end,2))-min(i_cndwt_1(i:end,2));

i=round(0.6/dt);
u_amplitude(2)=max(i_cndwt_2(i:end,2))-min(i_cndwt_2(i:end,2));

i=round(0.7/dt);
u_amplitude(3)=max(i_cndwt_3(i:end,2))-min(i_cndwt_3(i:end,2));

i=round(0.8/dt);
u_amplitude(4)=max(i_cndwt_4(i:end,2))-min(i_cndwt_4(i:end,2));

i=round(0.9/dt);
u_amplitude(5)=max(i_cndwt_5(i:end,2))-min(i_cndwt_5(i:end,2));

Tl=[0.5,0.6,0.7,0.8,0.9];

figure('position', [100, 100, 700, 225]);
plot(Tl,u_amplitude,'LineStyle',line_style{1},'Color',line_color);
grid on;
box off;

%% Task 4
% -------------------------------------------------------------

tde=csvread('task_dynamic_explicit.log');
tdep=csvread('task_dynamic_explicit_parallel.log');

figure('position', [100, 100, 700, 225]);
hold on;
plot(tde(:,1),tde(:,2),'LineStyle',line_style{1},'Color',line_color,'DisplayName','Explicit');
plot(tdep(:,1),tdep(:,2),'LineStyle',line_style{2},'Marker','o','Color',line_color,'DisplayName','Explicit Parallel');
grid on;
box off;
ylabel('u_{2}(m)');
xlabel('L (m)');
legend('show');

%% Task 5
% -------------------------------------------------------------

tdi=csvread('task_dynamic_implicit.log');
tdip=csvread('task_dynamic_implicit_parallel.log');

figure('position', [100, 100, 700, 225]);
hold on;
plot(tdi(:,1),tdi(:,2),'LineStyle',line_style{1},'Color',line_color,'DisplayName','Implicit');
plot(tdip(:,1),tdip(:,2),'LineStyle',line_style{2},'Marker','o','Color',line_color,'DisplayName','Implicit Parallel');
grid on;
box off;
ylabel('u_{2}(m)');
xlabel('L (m)');
legend('show');