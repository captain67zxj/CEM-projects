clear;clc;
%Extract the mesh data generated from pdetool:
load Microstrip_mesh.mat    %load mesh data generated from pdetool
co = p.';                   %coordinates of each node
N_by4 = t.';             
N = N_by4(:,(1:3));         %node number of each element
edge = e.';                 %Edge # node1 node2, including both Dirichlet and Neumann boundary nodes

%Coefficients for each each element:
a1 = zeros(size(N,1),1); a2 =zeros(size(N,1),1); a3= zeros(size(N,1),1);
b1 = zeros(size(N,1),1); b2 =zeros(size(N,1),1); b3= zeros(size(N,1),1);
c1 = zeros(size(N,1),1); c2 =zeros(size(N,1),1); c3= zeros(size(N,1),1);
area = zeros(size(N,1),1);  %area of each element

%% Assign coordinates and coefficients, and calculate area for each element
for i = 1:size(N,1)  %size(N,1)= # of elements 
    x1 = co(N(i,1),1); y1 = co(N(i,1),2);
    x2 = co(N(i,2),1); y2 = co(N(i,2),2);
    x3 = co(N(i,3),1); y3 = co(N(i,3),2);
    
    b1(i) = y2-y3; c1(i) = x3-x2;
    b2(i) = y3-y1; c2(i) = x1-x3;                           %Eq.(19)
    b3(i) = y1-y2; c3(i) = x2-x1;
    
    area(i) = 0.5*(b1(i)*c2(i)-b2(i)*c1(i));                %Eq.(20)
end
bb = [b1 b2 b3];
cc = [c1 c2 c3];
L = zeros(size(co,1),size(co,1));
epslion = 1;

%% Assmbly [L]
for i = 1:size(N,1)
    node = zeros(3,1);%local nodes
    node(1)=N(i,1); 
    node(2)=N(i,2); 
    node(3)=N(i,3);
       
     for l = 1:3 
         for k = 1:3                                      
              
             Int = epslion/(4*area(i))*(bb(i,l)*bb(i,k)+ cc(i,l)*cc(i,k)); %Intergral Eq.(22)

             L(node(l),node(k)) = L(node(l),node(k)) + Int;
                    
         end
     end
end

%% Assign the nodes on the dirichlet boundaries from e file and the mesh figure
      edge_D = [edge(1:30,1);edge(30,2);      % Dirichlet boundary with phi = 0
                edge(31:43,1);edge(43,2)];    % Dirichlet boundary with phi = 1000 V

%Assign the potential on the dirichlet boundaries
M = 32;
phi_D = zeros(size(edge_D,1),1);
for i = M:size(edge_D,1)
     phi_D(edge_D(i)) = 1000;       %High potential on the trace
end

%% Solve phi_unknown;
phi_unknown = - L(46:end,46:end)\(L(46:end,1:45)* phi_D);  % Solve phi_unknown
phi_all = [phi_D; phi_unknown];                            % Potential in the entire simulation domain 
U = scatteredInterpolant(co(:,1),co(:,2),phi_all);         % Interpolatation

%% triangular surface plot
trisurf(N,co(:,1),co(:,2),phi_all,'facecolor','interp')    %plot phi_all
zlabel('Potential (V)')
colorbar
ax = gca;
ax.FontSize = 24;
ax.LineWidth = 2;
ax.TickDir = 'both';
ax.TickLength = [0.008 0.01];
ax.FontWeight = 'bold';


%% plot the equipotential lines and static electric field distribution
% x = 0:0.01:0.8;
% y = 0:0.01:0.8;
% 
% phi_FEM_contour = zeros(length(x),length(y));
% for i = 1:length(x)
%     for j = 1:length(y)
%         phi_FEM_contour(i,j) = U(y(j),x(i));
%     end
% end
% 
% spacing = 0.01;
% [X,Y] = meshgrid(x,y);
% [DX,DY] = gradient(phi_FEM_contour,spacing);
% contourf(X,Y,phi_FEM_contour,10)
% colorbar
% hold on
% q = quiver(X,Y,-DX,-DY,'k','LineWidth',1);
% set(q,'AutoScale','on', 'AutoScaleFactor',2)
% xlabel('x (dm)')
% ylabel('y (dm)')
% axis equal
% ax = gca;
% ax.FontSize = 24;
% ax.LineWidth = 2;
% ax.TickDir = 'both';
% ax.TickLength = [0.008 0.01];
% ax.FontWeight = 'bold';
% hold off
