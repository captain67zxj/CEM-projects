clear;
clc;

%1D Lumped element transmission line Parameters:
C_0 = 0.298e-9;                             %[F/m] Capacitance
L_0 = 472e-9;                               %[H/m] Inductance
v = 1/(C_0*L_0)^0.5;                        %[m/s] Phase velocity
V_p = -10000;                               %[V] Input pulse voltage

%Ferrite parameters:
mu_0 = 4*pi*10^-7;                          %[H/m] Permeability of free space
d_i = 1.2e-3;                               %[m] Ferrite inner diameter
d_m = 3.5e-3;                               %[m] Ferrite outer diameter
alpha = 0.1;                                %Damping factor
gamma = 1.76e11;                            %[rad/s/T] Gyromagnetic ratio   
M_s = 0.35/mu_0;                            %[A/m] Saturation Magnetization
C_V = 0.5*mu_0*(d_m-d_i);                   %Constant in dV/dx
C_M = alpha*gamma*mu_0*M_s/(1+alpha^2);     %Constant in dM/dt
d_e = (d_m-d_i)/log(d_m/d_i);               %[m]Mean ferrite radius

%Simulation parameters:
L = 200e-3;                                 %[m] the tatal length of the simulation region/line length
dz = 200e-4;                                %[m] space step
n = fix(L/dz);                              %Number of sections
T = 10e-9;                                  %[s] Total simulation time
dt = 0.01*dz*(L_0*C_0)^0.5;                 %[s] Time step
C_C = dt/(dz*C_0);                          %Constant in Leapfrog
C_L = dt/(dz*L_0);                          %Constant in Leapfrog
z = linspace(0,L,n);

V = zeros(n,1);                             %voltage
Vp=V;                                       %V at time k+1
I = zeros(n,1);                             %Current
Ip=I;                                       %I at time k+1
M = zeros(n,1);                             %Magnetization
Mp=M;                                       %M at time k+1

V_input = zeros(fix(T/dt),1);               %Input pulse voltage varying with time
V_output = zeros(fix(T/dt),1);              %Onput pulse voltage varying with time
M_output= zeros(fix(T/dt),1);               %Magnetization varying with time

mu = 3e-9;                                  %Fitting parameter in the input voltage
a = 1e-10;                                  %Fitting parameter in the input voltage

%Initial conditions
t = 0;
tt = 0;
while (t<T)
    
    %Input voltage pulse, step function
    Vp(1) = V_p*(1-1/(exp((t-mu)/a)+1));                         
    
    %Input voltage pulse used for validate the linear line simulation:    
%     Vp(1) = V_p*sin(10*pi*t/T);                                            

    Vp(end) = V(end-1)+(dz-v*dt)/(dz+v*dt)*(V(end)-Vp(end-1));   %ABC

    t=t+dt;
    V = Vp; 
    I = Ip;
    M = Mp;

    for i = 2:n
        
        Mp(i) = dt*C_M*I(i)*(1-(M(i)/M_s)^2)/(pi*d_e)-M(i);
        
        Ip(i) = I(i)-C_L*(V(i)-V(i-1))-C_V*(Mp(i)-M(i))/L_0;
        
    end
    
    for i = 1:n-1
        
        Vp(i) = V(i)-C_C*(Ip(i+1)-Ip(i));

    end
    
    
%     %Voltage pulse papagating through space  
% 
%     clf;
%     plot(z,Vp);
%     title(sprintf('t= %.2f', t*1.e9));
% %     axis([0 L 0 -1.2*V_p])
%     xlabel('L')
%     ylabel('Voltage (V)')
%     ylim([-12000 12000])
%     shg; pause(0.01);
%     grid on;

%Save the voltage at the two ends (input and output) as a function of time
    tt= tt+1;
    V_input(tt) = V(1); 
    V_output(tt) = V(end);
    M_output(tt) = Mp(end);
end

% %##0 plot V and M
% subplot(1,2,1)
% m = 1:1:fix(T/dt);
% plot(m*dt,V_input(m),'LineWidth',2,'MarkerSize',12)
% hold on;
% plot(m*dt,V_output(m),'LineWidth',2,'MarkerSize',12)
% hold on;
% grid on; 
% legend({'Input - 1D FDTD','Output - 1D FDTD'},'box','off','Location','northeast')
% xlabel('Time (s)')
% ylabel('Voltage (V)')
% % ylim([0 1.2*V_p])
% xlim([0 T])
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.8]);
% ax = gca;
% ax.FontSize = 24;
% ax.LineWidth = 2;
% ax.TickDir = 'both';
% ax.TickLength = [0.008 0.01];
% ax.FontWeight = 'bold';
% 
% subplot(1,2,2)
% plot(m*dt,M_output(m),'LineWidth',2,'MarkerSize',12)
% grid on; 
% % legend({'Input - 1D FDTD','Output - 1D FDTD'},'box','off','Location','northeastoutside')
% xlabel('Time (s)')
% ylabel('M')
% % ylim([-50 1.2*V_p])
% % xlim([0 20e-9])
% 
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.8]);
% ax = gca;
% ax.FontSize = 24;
% ax.LineWidth = 2;
% ax.TickDir = 'both';
% ax.TickLength = [0.008 0.01];
% ax.FontWeight = 'bold';
% 
% 


% #1 Plot V without the data from COMSOL
m = 1:1:fix(T/dt);
plot(m*dt,V_input(m),'LineWidth',2,'MarkerSize',12)
hold on;
plot(m*dt,V_output(m),'LineWidth',2,'MarkerSize',12)
hold on;

grid on; 
legend({'Input - 1D FDTD','Output - 1D FDTD'},'box',...
    'off','Location','northeast')
xlabel('Time (s)')
ylabel('Voltage (V)')
xlim([0 T])

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.8]);
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 2;
ax.TickDir = 'both';
ax.TickLength = [0.008 0.01];
ax.FontWeight = 'bold';




% % #2 Plot V with the data from COMSOL
% filename = 'GNLTL_COMSOL.csv';
% num = xlsread(filename);
% t_c = num(:,1);
% V_in = num(:,2);
% V_out = num(:,3);
% N = 160;
% 
% m = 1:1:fix(T/dt);
% plot(m*dt,V_input(m),'LineWidth',2,'MarkerSize',12)
% hold on;
% plot(m*dt,V_output(m),'LineWidth',2,'MarkerSize',12)
% hold on;
% 
% plot(t_c(1:1:N)*10^-9, V_in(1:1:N),'-.','LineWidth',2,'MarkerSize',12)
% hold on
% plot(t_c(1:1:N)*10^-9, V_out(1:1:N),'-.','LineWidth',2,'MarkerSize',12)
% 
% grid on; 
% legend({'Input - 1D FDTD','Output - 1D FDTD','Input - 3D COMSOL',...
%     'Output - 3D COMSOL'},'box','off','Location','northeast')
% xlabel('Time (s)')
% ylabel('Voltage (V)')
% % ylim([0 1.2*V_p])
% xlim([0 T])
% 
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.8]);
% ax = gca;
% ax.FontSize = 24;
% ax.LineWidth = 2;
% ax.TickDir = 'both';
% ax.TickLength = [0.008 0.01];
% ax.FontWeight = 'bold';



