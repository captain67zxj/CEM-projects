clear;
clc;

%Simulation parameters:
c = 3e9;
L = 400e-3;                             %[m] the tatal length of the simulation region/line length
dx = 200e-5;                            %[m] space step
n = fix(L/dx);                          % Number of sections
T = 20e-9;                              %[s] Total simulation time
x = linspace(0,L,n);

C_0 = 1.5014e-10;                       %[F/m] Capacitance 
L_0 = 3.7379e-7;                        %[H/m] Inductance
V_p = 20000;

V = zeros(n,1);   %Voltage
Vp=V; %V at time n+1

I = zeros(n,1);   %Voltage
Ip=I; %V at time n+1

dt = dx*(L_0*C_0)^0.5;     

CFL_C = dt/(dx*C_0);
CFL_L = dt/(dx*L_0);

t = 0;
while (t<T)
    
%     Vp([1 end])=0;%Reflecting BC
      %Vp(1) = -(c*dt/dx)*(V(2)-V(1))+V(1);
     %Vp(end) = -(c*dt/dx)*(V(end)-V(end-1))+V(end-1);
   Vp(end) = V(end-1);
    Vp(1) = V_p*sin(40*pi*t/T);
      
    t=t+dt;
    V = Vp; 
    I = Ip;
    %save current and previous arrays (when t+dt, we assign Vm1 = V...
    %and assign Vp1 to V)
    
    %%source
%     V(1) = 1*gauspuls(t,1,0.1);

    for i = 2:n
        
        Ip(i) = I(i)+CFL_L*(V(i)-V(i-1));
        
    end

    for i = 2:n-1 
        
        Vp(i) = V(i)+CFL_C*(Ip(i+1)-Ip(i));
          
    end
    
    clf;
    %t
    %V(1)
    plot(x,Ip);
    title(sprintf('t= %.2f', t*1.e9));
%     axis([0 L -2*V_p 2*V_p])
    shg; pause(0.01);
end

%     plot(t,V(end-1));
    