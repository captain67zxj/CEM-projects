clear;clc;
%Extract the mesh data generated from pdetool:
load coax_mesh.mat
co = p.';          %coordinates of each node    .size(co,1) = # of nodes = 138
N_by4 = t.';             
N = N_by4(:,(1:3));%node number of each element. size(N,1)= # of elements = 226
edge = e.';        %dirichlet edge # node1 node2

%Coefficients for each each element:
a1 = zeros(size(N,1),1); a2 =zeros(size(N,1),1); a3= zeros(size(N,1),1);
b1 = zeros(size(N,1),1); b2 =zeros(size(N,1),1); b3= zeros(size(N,1),1);
c1 = zeros(size(N,1),1); c2 =zeros(size(N,1),1); c3= zeros(size(N,1),1);
%area of each element
area = zeros(size(N,1),1); 

%Assign coordinates and coefficients, and calculate area for each element
for i = 1:size(N,1)  %size(N,1)= # of elements 
    x1 = co(N(i,1),1); y1 = co(N(i,1),2);
    x2 = co(N(i,2),1); y2 = co(N(i,2),2);
    x3 = co(N(i,3),1); y3 = co(N(i,3),2);
    
    a1(i) = x2*y3-x3*y2; b1(i) = y2-y3; c1(i) = x3-x2;
    a2(i) = x3*y1-x1*y3; b2(i) = y3-y1; c2(i) = x1-x3;      %Eq.(20)
    a3(i) = x1*y2-x2*y1; b3(i) = y1-y2; c3(i) = x2-x1;
    
    area(i) = 0.5*(b1(i)*c2(i)-b2(i)*c1(i));                %Eq.(21)
end

bb = [b1 b2 b3];
cc = [c1 c2 c3];

phi_unknown = zeros(size(co,1)-size(edge,1),1);                 %unknown potential
K = zeros(size(co,1),size(co,1));
L = zeros(size(co,1)-size(edge,1),size(co,1)-size(edge,1));
% b = zeros(size(co,1)-size(edge,1),1);
epslion = 1;%8.854e-12;

for i = 1:size(N,1)
    node = zeros(3,1);%local nodes
    node(1)=N(i,1); node(2)=N(i,2); node(3)=N(i,3);
    
    
     for l = 1:3 
         for k = 1:3
             
                          
              if  any(edge(:,1)  == node(k)) || any(edge(:,1) == node(l))
                  continue;
              end
              
             Int = epslion/(4*area(i))*(bb(i,l)*bb(i,k)+ cc(i,l)*cc(i,k)); %intergral
             if (Int == 0)
                 disp(Int);
             end
             %disp(Int);
             K(node(l),node(k)) = K(node(l),node(k)) + Int;
              
             
         end
     end
end


ii = 1;
jj = 1;
for i = 1: size(co,1)
    jj = 1;
    if  any(edge(:,1) == i)
        continue;
    end
    
    for j = 1:size(co,1)
        if  any(edge(:,1) == j)
            continue;
        end
        L(ii,jj) = K(i,j);
        jj = jj +1;
    end
    
    ii = ii + 1;
end
nnz(K)
issymmetric(K)
nnz(L)
issymmetric(L)

% phi = ones(size(KK,1),1);
M = 30;
phi_D = zeros(size(edge,1),1);
for i = M:size(edge,1)
     phi_D(edge(i)) = 1000;       %High potential on the inner boundary
end

phi_unknown = - L(51:end,51:end)\(L(51:end,1:50)* phi_D);  % Solve phi_unknown
phi_all = [phi_D; phi_unknown];                            % Potential in the entire simulation domain 

trisurf(N,co(:,1),co(:,2),phi_all,'facecolor','interp')    %plot phi_all

