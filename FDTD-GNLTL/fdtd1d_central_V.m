clear;
clc;

Lx = 10;
dx = 0.1;
nx = fix(Lx/dx);

x = linspace(0,Lx,nx);

T = 50;

V = zeros(nx,1);   %Voltage
Vm=V; %V at time n-1
Vp=V; %V at time n+1

% I = zeros(nx,1);   %Current
% Im=I; %V at time n-1
% Ip=V; %V at time n+1

CFL = 1;
% b=0.01;
b=0;
c=1;

dt=CFL*dx/c;


t = 0;
while (t<T)
    
    
%     V([1 end])=0; %Reflecting BC
    Vp(1) = V(2)+((CFL-1)/(CFL+1))*(Vp(2)-V(1));%%Absorbing BC
    Vp(end) = V(end-1)+((CFL-1)/(CFL+1))*(Vp(end-1)-V(end));%%Absorbing BC

    %solution
    t=t+dt;
    Vm = V; V = Vp; 
    %save current and previous arrays (when t+dt, we assign Vm1 = V...
    %and assign Vp1 to V)
    
    %%source
    V(1) = dt^2*20*sin(20*pi*t/T);
%     V(1) = 1*gauspuls(t,1,0.1);

    for i = 2:nx-1
        Vp(i) = 2*V(i)- Vm(i) + CFL^2*(V(i+1)-2*V(i)+V(i-1));
    end
    
    
    clf;
    plot(x,V);
    title(sprintf('t= %.2f', t));
    axis([0 Lx -0.5 0.5])
    shg; pause(0.01);
end

%     plot(t,V(end-1));
    