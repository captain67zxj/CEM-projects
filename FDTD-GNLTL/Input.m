clear;
clc;

%1D Lumped element transmission line Parameters:
mu_0 = 4*pi*10^-7;                            %[H/m] Permeability of free space
c = 3e9;                                    %Speed of light
C_0 = 0.298e-9;                             %[F/m] Capacitance
L_0 = 472e-9;                            %[H/m] Inductance
V_p = -10000;                                %[V] Input pulse voltage

%Ferrite parameters:
d_i = 1.2e-3;%3e-3;                                 %[m] Ferrite inner diameter
d_m = 3.5e-3;%5e-3;                                 %[m] Ferrite outer diameter
alpha = 0.1;                                %Damping factor
gamma = 1.76e11;                            %[rad/s/T] Gyromagnetic ratio   
M_s = 0.35/mu_0;                                  %[A/m] Saturation Magnetization
C_V = 0.5*mu_0*(d_m-d_i);                     %Constant in dV/dx
C_M = alpha*gamma*mu_0*M_s/(1+alpha^2);     %Constant in dM/dt
d_e = (d_m-d_i)/log(d_m/d_i);               %[m]Mean ferrite radius

%Simulation parameters:
L = 200e-3;                                 %[m] the tatal length of the simulation region/line length
dz = 200e-5;                                %[m] space step
n = fix(L/dz);                              %Number of sections
T = 20e-9;                                  %[s] Total simulation time

%Other parameters:
dt = 0.1*dz*(L_0*C_0)^0.5;                      %[s] time step
C_C = dt/(dz*C_0);                          %Constant
C_L = dt/(dz*L_0);                          %Constant

z = linspace(0,L,n);

V = zeros(n,1);                             %voltage
Vp=V;
%V at time k+1
I = zeros(n,1);                             %Current
Ip=I;
%I at time k+1
M = zeros(n,1);                             %Magnetization
Mp=M;                                       %M at time k+1

V_input = zeros(fix(T/dt),1);               %Input pulse voltage varying with time
V_output = zeros(fix(T/dt),1);              %Onput pulse voltage varying with time

M_output= zeros(fix(T/dt),1);

mu = 3e-9;
a = 1e-10;

%Initial conditions
t = 0;
tt = 0;
while (t<T)
    
%      Vp(1) = -(c*dt/dz)*(V(2)-V(1))+V(1);
%     Vp(end) = 0;

    Vp(end) = V(end-1);                     %ABC
    Vp(1) = V_p*(1-1/(exp((t-mu)/a)+1));    %Input pulse
%     Mp(1) = 25000;

    t=t+dt;
    V = Vp; 
    I = Ip;
    M = Mp;

    
    for i = 2:n
        
        Mp(i) = dt*C_M*I(i)*(1-(M(i)/M_s)^2)/(pi*d_e)-M(i);
              
        Ip(i) = I(i)-C_L*(V(i)-V(i-1))-C_V/L_0*(Mp(i)-M(i));
        
    end
    
    for i = 1:n-1
        
        Vp(i) = V(i)-C_C*(Ip(i+1)-Ip(i));
          
    end
    
%        Mp(end)
    
%     %Voltage pulse papagating through space  

%     clf;
%     plot(z,Vp);
%     title(sprintf('t= %.2f', t*1.e10));
% %     axis([0 L 0 1.2*V_p])
%     shg; pause(0.01);    


    %Save the voltage at the two ends (input and output) as a function of time
    tt= tt+1;
    V_input(tt) = V(1); 
    V_output(tt) = V(end-1);
    M_output(tt) = Mp(50);

end



%Plot input and output voltage as a function of time
% subplot(1,2,1)

m = 1:1:fix(T/dt);
plot(m*dt,V_input(m),'LineWidth',2,'MarkerSize',12)
hold on;
plot(m*dt,V_output(m),'LineWidth',2,'MarkerSize',12)
hold on;


grid on; 
legend({'Input - 1D FDTD','Output - 1D FDTD'},'box','off','Location','northeastoutside')
xlabel('Time (s)')
ylabel('Voltage (V)')
% ylim([0 1.2*V_p])
xlim([0 20e-9])

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.8]);
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 2;
ax.TickDir = 'both';
ax.TickLength = [0.008 0.01];
ax.FontWeight = 'bold';




% subplot(1,2,2)
% plot(m*dt,M_output(m),'LineWidth',2,'MarkerSize',12)
% grid on; 
% % legend({'Input - 1D FDTD','Output - 1D FDTD'},'box','off','Location','northeastoutside')
% xlabel('Time (s)')
% ylabel('M')
% % ylim([-50 1.2*V_p])
% % xlim([0 20e-9])
% 
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 0.8]);
% ax = gca;
% ax.FontSize = 24;
% ax.LineWidth = 2;
% ax.TickDir = 'both';
% ax.TickLength = [0.008 0.01];
% ax.FontWeight = 'bold';