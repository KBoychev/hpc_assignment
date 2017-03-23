import math
import matplotlib.pyplot as plt
import numpy as np
import csv


def load_csv_file(csv_file_name = ""):   

   csv_file_contents_1=[]
   csv_file_contents_2=[]

   with open(csv_file_name, "rb") as csv_file: 	
	csv_file_reader = csv.reader(csv_file)
	for csv_row in csv_file_reader:
		csv_file_contents_1.append(float(csv_row[0]))
		csv_file_contents_2.append(float(csv_row[1]))
   return csv_file_contents_1, csv_file_contents_2;

## Task1
# -------------------------------------------------

E=210000000000
I=0.0000144
L=10
qy=1000
Fy=1000

x,u2=load_csv_file("task_static.log")

u2_a=[]
for i in range(0,len(x)):
	u2_a.append((Fy*pow(x[i],2))/(48*E*I)*(3*L-4*x[i])+(qy*pow(x[i],2))/(24*E*I)*pow(L-x[i],2))


plt.figure(0)
plt.plot(x,u2,linestyle='-',color = '#1c4587',label='Numerical')
plt.plot(x,u2_a,linestyle='none',marker='+',color = '#1c4587',label='Analytical')
plt.xlabel('x (m)')
plt.ylabel('u2 (m)')
plt.legend()
plt.grid(True)
plt.savefig("task1.png")



## Task2
# -------------------------------------------------



t,u2_T=load_csv_file("explicit_center_node_defletion_wtr_time_TleqT.log")
t,u2_05T=load_csv_file("explicit_center_node_defletion_wtr_time_Tleq05T.log")
t,u2_06T=load_csv_file("explicit_center_node_defletion_wtr_time_Tleq06T.log")
t,u2_07T=load_csv_file("explicit_center_node_defletion_wtr_time_Tleq07T.log")
t,u2_08T=load_csv_file("explicit_center_node_defletion_wtr_time_Tleq08T.log")
t,u2_09T=load_csv_file("explicit_center_node_defletion_wtr_time_Tleq09T.log")

plt.figure(1)
plt.plot(t,u2_T,linestyle='-',color = '#1c4587',label='Tl=T')
plt.plot(t,u2_05T,linestyle='--',color = '#1c4587',label='Tl=0.5T')
plt.plot(t,u2_06T,linestyle='-.',color = '#1c4587',label='Tl=0.6T')
plt.plot(t,u2_07T,linestyle=':',color = '#1c4587',label='Tl=0.7T')
plt.plot(t,u2_08T,linestyle='none',marker='+',markevery=80,color = '#1c4587',label='Tl=0.8T')
plt.plot(t,u2_09T,linestyle='none',marker='x',markevery=80,color = '#1c4587',label='Tl=0.9T')
plt.xlabel('t (s)')
plt.ylabel('u2 at L/2 (m)')
plt.legend()
plt.grid(True)
plt.savefig("task2.png")



T=1.0
Nt=8000.0
dt=T/(Nt-1.0)

u2_amp=[]


i=int(round(0.5/dt))
u2_amp.append(max(u2_05T[i:])-min(u2_05T[i:]))

i=int(round(0.6/dt))
u2_amp.append(max(u2_06T[i:])-min(u2_06T[i:]))

i=int(round(0.7/dt))
u2_amp.append(max(u2_07T[i:])-min(u2_07T[i:]))

i=int(round(0.8/dt))
u2_amp.append(max(u2_08T[i:])-min(u2_08T[i:]))

i=int(round(0.9/dt))
u2_amp.append(max(u2_09T[i:])-min(u2_09T[i:]))

Tl=[0.5,0.6,0.7,0.8,0.9];

plt.figure(2)
plt.loglog(Tl,u2_amp,linestyle='-',color = '#1c4587')
plt.xlabel('Tl (s)')
plt.ylabel('amplitude of u2 at L/2 (m)')
plt.grid(True)
plt.savefig("task2_amplitude.png")


## Task3
# -------------------------------------------------


t,u2_T=load_csv_file("implicit_center_node_defletion_wtr_time_TleqT.log")
t,u2_05T=load_csv_file("implicit_center_node_defletion_wtr_time_Tleq05T.log")
t,u2_06T=load_csv_file("implicit_center_node_defletion_wtr_time_Tleq06T.log")
t,u2_07T=load_csv_file("implicit_center_node_defletion_wtr_time_Tleq07T.log")
t,u2_08T=load_csv_file("implicit_center_node_defletion_wtr_time_Tleq08T.log")
t,u2_09T=load_csv_file("implicit_center_node_defletion_wtr_time_Tleq09T.log")

plt.figure(3)
plt.plot(t,u2_T,linestyle='-',color = '#1c4587',label='Tl=T')
plt.plot(t,u2_05T,linestyle='--',color = '#1c4587',label='Tl=0.5T')
plt.plot(t,u2_06T,linestyle='-.',color = '#1c4587',label='Tl=0.6T')
plt.plot(t,u2_07T,linestyle=':',color = '#1c4587',label='Tl=0.7T')
plt.plot(t,u2_08T,linestyle='none',marker='+',markevery=10,color = '#1c4587',label='Tl=0.8T')
plt.plot(t,u2_09T,linestyle='none',marker='x',markevery=10,color = '#1c4587',label='Tl=0.9T')
plt.xlabel('t (s)')
plt.ylabel('u2 at L/2 (m)')
plt.legend()
plt.grid(True)
plt.savefig("task3.png")



T=1.0
Nt=1000.0
dt=T/(Nt-1.0)

u2_amp=[]


i=int(round(0.5/dt))
u2_amp.append(max(u2_05T[i:])-min(u2_05T[i:]))

i=int(round(0.6/dt))
u2_amp.append(max(u2_06T[i:])-min(u2_06T[i:]))

i=int(round(0.7/dt))
u2_amp.append(max(u2_07T[i:])-min(u2_07T[i:]))

i=int(round(0.8/dt))
u2_amp.append(max(u2_08T[i:])-min(u2_08T[i:]))

i=int(round(0.9/dt))
u2_amp.append(max(u2_09T[i:])-min(u2_09T[i:]))

Tl=[0.5,0.6,0.7,0.8,0.9];

plt.figure(4)
plt.plot(Tl,u2_amp,linestyle='-',color = '#1c4587')
plt.xlabel('Tl (s)')
plt.ylabel('amplitude of u2 at L/2  (m)')
plt.grid(True)
plt.savefig("task3_amplitude.png")



## Task4
# -------------------------------------------------


x,u2=load_csv_file("task_dynamic_explicit.log")
x1,u2_p=load_csv_file("task_dynamic_explicit_parallel.log")

plt.figure(6)
plt.plot(x,u2,linestyle='-',color = '#1c4587',label='Serial')
plt.plot(x1,u2_p,linestyle='none',marker='+',color = '#1c4587',label='Parallel')
plt.xlabel('x (m)')
plt.ylabel('u2 (m)')
plt.legend()
plt.grid(True)
plt.savefig("task4.png")



## Task5
# -------------------------------------------------

x,u2=load_csv_file("task_dynamic_implicit.log")
x,u2_p=load_csv_file("task_dynamic_implicit_parallel.log")

plt.figure(7)
plt.plot(x,u2,linestyle='-',color = '#1c4587',label='Serial')
plt.plot(x,u2_p,linestyle='none',marker='+',color = '#1c4587',label='Parallel')
plt.xlabel('x (m)')
plt.ylabel('u2 (m)')
plt.legend()
plt.grid(True)
plt.savefig("task5.png")


