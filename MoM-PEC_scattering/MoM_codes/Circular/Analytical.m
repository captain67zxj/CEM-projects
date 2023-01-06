clear;clc;
j = sqrt(-1);
R = 0.5;                        %Radius of the circular cylinder
lambda = R;                     %Wavelength
% lambda = R/2;                   %Wavelength
k_0 = 2*pi/lambda;
%% Calculating incident field, scattering field, and total field
x = linspace(-5*R,5*R); y = linspace(-5*R,5*R); 
phi = zeros(size(x,2),size(y,2));
Inc_field = zeros(size(x,2),size(y,2));
E_total = zeros(size(x,2),size(y,2)); E_sc = zeros(size(x,2),size(y,2));
H_total = zeros(size(x,2),size(y,2)); H_sc = zeros(size(x,2),size(y,2));
N = 80; 
[X,Y] = meshgrid(x,y);

for i = 1:size(x,2)
    
    for ii = 1:size(y,2)
        if x(i)^2+y(ii)^2 < R^2
            H_sc(i,ii) = 0 ; E_total (i,ii) = 0;
            E_sc(i,ii) = 0 ; H_total (i,ii) = 0;
        else
            
        for p = -N:N
            phi(i,ii) = atan(y(ii)/x(i));
            
            if x(i) < 0
                phi(i,ii) = phi(i,ii) + pi;
            end
            
            Inc_field(i,ii) = Inc_field(i,ii) + real(j^-p*besselj(p,k_0*...
                sqrt(x(i)^2+y(ii)^2))*exp(j*p*phi(i,ii))); % Eq.(29)&(31)
            
            E_sc(i,ii) = E_sc(i,ii) - real(j^-p*besselj(p,k_0*R)*...
                besselh(p,2,k_0*sqrt(x(i)^2+y(ii)^2))*exp(j*p*phi(i,ii))/...
                besselh(p,2,k_0*R)) ;% Eq.(30)
            
            H_sc(i,ii) = H_sc(i,ii) - real(j^-p*(-besselj(p+1,k_0*R)+p/(k_0*R)*...
                besselj(p,k_0*R))*besselh(p,2,k_0*sqrt(x(i)^2+y(ii)^2))*exp(j*p*phi(i,ii))/...
                (-besselh(p+1,2,k_0*R)+p/(k_0*R)*besselh(p,2,k_0*R)));% Eq.(32)
        end
        E_total (i,ii) = Inc_field(i,ii) + E_sc(i,ii);
        H_total (i,ii) = Inc_field(i,ii) + H_sc(i,ii);
        end
    end

end

%% Plotting
%% Incident field
% contourf(X,Y,Inc_field',20)
% colorbar
% hold on
% axis equal
% xlabel('x ')
% ylabel('y ')
% ax = gca;
% ax.FontSize = 24;
% ax.TickDir = 'both';
% ax.TickLength = [0.008 0.01];
% ax.FontWeight = 'bold';
% hold off

%% scattering field and total field

p = nsidedpoly(1000, 'Center', [0 0], 'Radius', R);

subplot(2,2,1)
contourf(X,Y,E_sc',20)
colorbar
hold on
plot(p, 'FaceColor', 'w')
hold on
axis equal
xlabel('x ')
ylabel('y ')
title('E_z^{sc}')
ax = gca;
ax.TickDir = 'both';
ax.TickLength = [0.008 0.01];
ax.FontWeight = 'bold';
hold off

subplot(2,2,2)
contourf(X,Y,H_sc',20)
colorbar
hold on
plot(p, 'FaceColor', 'w')
hold on
axis equal
xlabel('x ')
ylabel('y ')
title('H_z^{sc}')
ax = gca;
ax.TickDir = 'both';
ax.TickLength = [0.008 0.01];
ax.FontWeight = 'bold';
hold off

subplot(2,2,3)
contourf(X,Y,E_total',20)
colorbar
hold on
plot(p, 'FaceColor', 'w')
hold on
axis equal
xlabel('x ')
ylabel('y ')
title('E_z')
ax = gca;
ax.TickDir = 'both';
ax.TickLength = [0.008 0.01];
ax.FontWeight = 'bold';
hold off

subplot(2,2,4)
contourf(X,Y,H_total',20)
colorbar
hold on
plot(p, 'FaceColor', 'w')
hold on
axis equal
xlabel('x ')
ylabel('y ')
title('H_z')
ax = gca;
ax.TickDir = 'both';
ax.TickLength = [0.008 0.01];
ax.FontWeight = 'bold';
hold off
