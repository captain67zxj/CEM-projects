clear;clc;
%Extract the mesh data generated from pdetool:
load square_mesh.mat          %load mesh data generated from pdetool
co = p.';                   %coordinates of each node    .size(co,1) = # of nodes
edge = e.';                 %Dirichlet edge # node1 node2
%% extract and aligne the boundary edges into B_edge, each edge index
%associates two edge nodes
B_edge = zeros(size(edge,1),2);
for i = 1:size(edge,1)
    B_edge(i,1) = edge(i,1); 
    B_edge(i,2) = edge(i,2);
end
M = size(B_edge,1);  % number of segments

%% nth segment length calculated from the coordinates of the boundary nodes
s = zeros(M,1);
for i = 1:M
   s(i) = sqrt((co(B_edge(i,1),1)-co(B_edge(i,2),1))^2+...
            ((co(B_edge(i,1),2)-co(B_edge(i,2),2))^2));
end
%% center of the mth (or nth) segment (r_m or r_n)
r = zeros(M,2);
for i = 1:size(r)
    r(i,1) = (co(B_edge(i,1),1)+co(B_edge(i,2),1))/2;
    r(i,2) = (co(B_edge(i,1),2)+co(B_edge(i,2),2))/2;
end

%% %%%Matrix assembly
Zmn_TM = zeros(M,M);
Zmn_TE = zeros(M,M);
e_const = 2.7183;
gamma = 1.781;
j = sqrt(-1);
D = 0.5;                        %Side length
lambda = D/2;                     %Wavelength
k_0 = 2*pi/lambda;
Z_0 = 377;                      %Impedance of free space
%%  Assembly Zmn_TM for the TM case Eq.(14)
for m = 1:M
    for n = 1:M
        x_m = r(m,1); x_n = r(n,1);
        y_m = r(m,2); y_n = r(n,2);        
        Zmn_TM(m,n) = 1/4*k_0*Z_0*s(n)*...
            besselh(0,2,(k_0*sqrt((x_m-x_n)^2+(y_m-y_n)^2)));
    end
        Zmn_TM(m,m) = 1/4*k_0*Z_0*s(m)*(1-j*(2/pi)*...
            log(k_0*gamma*s(m)/(4*e_const)));
end

%%  Assembly Zmn_TE for the TE case Eq.(24)
for m = 1:M
    for n = 1:M
        x_m = r(m,1); x_n = r(n,1); y_m = r(m,2); y_n = r(n,2);       
        Zmn_TE(m,n) = k_0*s(n)/(4*j)*...
            besselh(1,2,(k_0*sqrt((x_m-x_n)^2+(y_m-y_n)^2)))*...
            ((x_n*(x_m-x_n)+y_n*(y_m-y_n))/sqrt(x_n^2+y_n^2))/...
            sqrt((x_m-x_n)^2+(y_m-y_n)^2);
    end
        Zmn_TE(m,m) = -1/2;
end

%% incident field
V = zeros(M,1);
N = 40;
for m = 1:M
    x_m = r(m,1);
    V(m) = exp(-j*k_0*x_m);
end 

%% sovling the surface current density
J_zn = Zmn_TM \ V; %TM
J_st = Zmn_TE \ V; %TE

%% Calculating incident field, scattering field, and total field
x = linspace(-5*D/2,5*D/2); y = linspace(-5*D/2,5*D/2); 
Inc_field = zeros(size(x,2),size(y,2));
E_total = zeros(size(x,2),size(y,2)); E_sc = zeros(size(x,2),size(y,2));
H_total = zeros(size(x,2),size(y,2)); H_sc = zeros(size(x,2),size(y,2));

[X,Y] = meshgrid(x,y);

for i = 1:size(x,2)
    Inc_field(i,:) = cos(-k_0*x(i));

    for j = 1:size(y,2) 
        if  -0.25 < x(i) && x(i) < 0.25 && -0.25 < y(j) && y(j)< 0.25
            H_sc(i,j) = 0; E_sc(i,j) = 0 ;
            E_total(i,j) = 0; H_total(i,j) = 0;
        else
            for m = 1:M
                x_m = r(m,1); y_m = r(m,2);
             
                %TM:
                G_0 = 1/(4*j)*besselh(0,2,(k_0*sqrt((x(i)-x_m)^2+(y(j)-y_m)^2)));
                E_sc(i,j) =  E_sc(i,j) - real(j*k_0*Z_0*G_0*J_zn(m)*s(m)) ;
             
                %TE:
                G_0_n = k_0/(4*j)*besselh(1,2,(k_0*sqrt((x(i)-x_m)^2+(y(j)-y_m)^2)))*...
                    ((x_m*(x(i)-x_m)+y_m*(y(j)-y_m))/sqrt(x_m^2+y_m^2))/...
                    sqrt((x(i)-x_m)^2+(y(j)-y_m)^2);
                H_sc(i,j) =  H_sc(i,j) - real(j*G_0_n*J_st(m)*s(m));
            end
            E_total (i,j) = cos(-k_0*x(i)) + E_sc(i,j);
            H_total (i,j) = cos(-k_0*x(i)) + H_sc(i,j); 
        end   
    end
end

%% Plotting
%%% Incident field
% contourf(X,Y,Inc_field',20)
% colorbar
% hold on
% axis equal
% xlabel('x ')
% ylabel('y ')
% ax = gca;
% ax.FontSize = 24;
% % ax.LineWidth = 2;
% ax.TickDir = 'both';
% ax.TickLength = [0.008 0.01];
% ax.FontWeight = 'bold';
% hold off

%% scattering field and total field
p = polyshape([-D/2 -D/2 D/2 D/2],[D/2 -D/2 -D/2 D/2]);

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

