% Rachel Murphy
%1351620
% Hemholtz Equation-Gauss Seidel Method
clc
clear
N=10; % The number of increments for both x and y dimensions
h=(2*pi)/(N-1); % The increment of both x and y dimensions
lambda=-pi; % Given constant

% Rectangle lengths
a_x=-pi;
a_y=a_x;
b_x=pi;
b_y=b_x;

x=a_x:h:b_x; % x distances
y=a_y:h:b_y; % y distances

% Creating the F matrix
F=zeros(N,N);
for i=1:N
    for j=1:N
        F(i,j)=cos((pi/2)*(2*((x(j)-a_x)/(b_x-a_x))+1))*sin(pi*(y(i)-a_y)/(b_y-a_y));
    end
end
F=lambda*h*F;

% Creating initial U_old matrix, which contains initial solution guesses
U_old=zeros(N,N);
U_old(1:N,1)=cos (pi*(y-a_y)-1).*cosh(b_y-y);
U_old(1:N,N)=(y-a_y).^2.*sin(pi*(y-a_y)/(2*(b_y-a_y)));

% Checkpointing

% Before the start of the iteration loop, "check-in" each variable
% that should be checkpointed in the event of restarting the job

matfile = 'HELM_GAUSS_SEIDEL.mat';     % mandatory; name of checkpoint mat-file
s = struct();                                % mandatory; create struct for checkpointing
s = chkin(s,{'iter'});                       % mandatory; iter is iteration loop index
s = chkin(s,{'frequency'});                  % mandatory; frequency is checkpointing period 
                                             % i.e., how often to perform a save


chkNames = fieldnames(s);    % the full list of variables to checkpoint
nNames = length(chkNames);   % number of variables in list

tic; % Timer to find the run time

for i=1:10
    % Iterative process to solve for the new matrix using the solutions from
    U_new=U_old;
    
    % Solves for boundary conditions of the lower edge using ghost
    % nodes/neumann conditions
    for i=2:N-1
        U_new(1,i)=(U_new(1,i+1)+U_new(1,i-1)+U_new(2,i)+U_new(2,i)-h^2*F(1,i))/(4-h^2*lambda);
    end
    
    % Solves for boundary conditions of the upper edge using ghost
    % nodes/neumann conditions
    
    for i=2:N-1
        U_new(N,i)=(U_new(N,i+1)+U_new(N,i-1)+U_new(N-1,i)+U_new(N-1,i)-h^2*F(1,i))/(4-h^2*lambda);
    end
    
    % Solves for the internal nodes (the numbers in between the boundaries)
    for i=2:N-1
        for j=2:N-1
            U_new(i,j)=(U_new(i+1,j)+ U_new(i-1,j)+U_new(i,j+1)+U_new(i,j-1)-h^2*F(i,j))/(4-h^2*lambda);
        end
    end
   U_error=abs((mean(mean(U_new))-mean(mean(U_old)))./(mean(mean(U_new))))*100 
   U_old=U_new;
end

runtime=toc; % End of timer

FinalValue = mean(mean(U_new (1:N,1:N))) % Needed for Grid Convergence Study

% Plotting the solution
figure(1)
mesh(x,y,U_new)                    % Surface Plot
xlabel('X','fontSize',12);
ylabel('Y','fontSize',12);
zlabel('U (x,y)')
title('Gauss-Seidel Method with 100 Nodes','fontsize',12);

figure(2)
contourf(x,y,U_new),colorbar 
xlabel('X','fontSize',12);
ylabel('Y','fontSize',12);
title('Contour Plot with Gauss-Seidel Method','fontsize',12)