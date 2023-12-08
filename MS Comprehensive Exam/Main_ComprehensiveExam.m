%% Clear
clear; close all; clc

%% Problem Givens
N = 50; %nodes
dt = .0005;
totalTime = .05; %in seconds

%rod params
rodLength = .1;
deltaL = rodLength/(N - 1);
b = .01;
h = .002;
Young = 200e9;
rho = 7850; %this only affects the time it takes to reach steady-state.

PLoad = [10]; %point load, in N
%PLoad = (0:5:100)'; %loop over various P
PSteps = size(PLoad);

%utility vars
ne = N - 1; %num of edges
I = b*(h^3)/12;
A = b*h;
EI = Young*I; %bending stiffness, CHANGE TO RECTANGLE
EA = Young*A; %axial stiffness, CHANGE TO RECTANGLE

%% Vector Setup
%initial position matrix, each row is initial x and y of a node
nodes = zeros(N, 2);
for k = 1:N
    nodes(k,1) = (k-1)*deltaL; %initial x of node at x = 0
    nodes(k,2) = 0; %initial y of node, all = 0
end

%initial DOF vector
q0 = zeros(2*N, 1);
for k = 1:N
    q0(2*k - 1) = nodes(k,1); %x-coord
    q0(2*k) = nodes(k,2); %y-coord
end

%mass matrix M
M = zeros(2*N, 2*N);
for k = 1:N
    M(2*k-1, 2*k-1) = (b*h*rodLength*rho)/(N-1); %assigns on every other entry of diagonal
    M(2*k, 2*k) = M(2*k-1,2*k-1); %mass is repeated in M
end

%set new position and velocity
q = q0;
u = (q - q0)/dt;

tol = (EI/(rodLength^2))*1e-3;
nSteps = round(totalTime/dt);

%vectors to store data over time
yPosMax = zeros(nSteps, 1);
yPos(1) = q(2);

%point force vector, P
P = zeros(2*N,1);
PNode = N; % point load is applied on last node

%% Time Loop

yPosLoad = zeros(PSteps(1),1);

for p = 1:PSteps %time loop for P
    
    P(2*PNode) = -PLoad(p);
    
    for c = 1:nSteps
        
        fprintf('Time = %f\n', (c-1)*dt);
        
        %boundary conditions
        fixedDOF = [1; 2; 4]; % left end is fixed, first rotation = 0 so y2 = 0
        freeDOF = [3; (5:(2*N))']; % all other nodes are free
        fixedBC = [0; 0; 0];
        
        q = q0; %initial guess
        q(fixedDOF) = fixedBC;
        
        %---Newton-Raphson loop---
        err = 10*tol;
        
        while (err > tol)
            
            qFree = q(freeDOF);
            
            %Calculating f and J....
            %inertial term
            f = M / dt * ( (q-q0) / dt - u );
            J = M / dt^2;
            
            %elastic terms
            %linear spring forces
            for k = 1:(N-1)
                xk = q(2*k - 1);
                yk = q(2*k);
                xkp1 = q(2*k + 1);
                ykp1 = q(2*k + 2);
                dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
                dJ = hessEs(xk, yk, xkp1, ykp1, deltaL, EA);
                ind = [2*k - 1; 2*k; 2*k + 1; 2*k + 2]; %sets to correct 4x4 spot in matrix
                f(ind) = f(ind) + dF;
                J(ind,ind) = J(ind,ind) + dJ;
            end
            
            % Bending spring forces
            for k = 2:(N-1)
                xkm1 = q(2*k - 3);
                ykm1 = q(2*k - 2);
                xk = q(2*k - 1);
                yk = q(2*k);
                xkp1 = q(2*k + 1);
                ykp1 = q(2*k + 2);
                curvature0 = 0;
                dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
                dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
                ind = [2*k - 3; 2*k - 2; 2*k - 1; 2*k; 2*k + 1; 2*k + 2];
                f(ind) = f(ind) + dF;
                J(ind,ind) = J(ind,ind) + dJ;
            end
            
            %point load term
            f = f - P;
            
            %Newton-Raphson calculation
            fFree = f(freeDOF);
            JFree = J(freeDOF, freeDOF);
            
            dqFree = JFree\fFree;
            qFree = qFree - dqFree;
            
            err = sum(abs(fFree));
            
            q(freeDOF) = qFree;
            
        end %----end of Newton-Raphson loop----
        
        %update old values to new
        u = (q - q0)/dt;
        q0 = q;
        
        %frame-by-frame plotting
        figure(1);
        plot( q(1:2:end), q(2:2:end), 'ro-');
        xlabel('x (m)')
        ylabel('y (m)')
        axis equal
        drawnow
        
        %store data
        yPos = zeros(N, 1);
        for k = 1:N
            yPos(k) = q(2*k);
        end
        yPosMax(c) = max(abs(yPos));
         
    end
    
    yPosLoad(p) = yPosMax(nSteps);
    
end

%% Plotting

timeArr = (1:nSteps)*dt;

% figure(2)
% hold on
% title('Absolute Maximum Vertical Displacement vs. Time')
% plot(timeArr, yPosMax)
% xlabel('time (s)');
% ylabel('maximum y-displacement (m)');
% hold off

EBdef = PLoad.*(rodLength^3)./(3*Young*I); %analytical Euler-Bernoulli beam deflection

% figure(2)
% hold on
% title('Absolute Maximum Vertical Displacement vs. Point Load')
% plot(PLoad, yPosLoad, 'DisplayName', 'Implicit Beam Simulation')
% plot(PLoad, EBdef, 'DisplayName', 'Euler-Bernoulli Beam Deflection')
% xlabel('point load (N)');
% ylabel('maximum y-displacement (m)');
% legend()
% hold off

 deviation = 100*abs((yPosLoad - EBdef)./EBdef);
% figure(3)
% hold on
% title('Deviation of Implicit Beam Simulation from Euler-Bernoulli Beam')
% plot(PLoad, deviation)
% xlabel('point load (N)');
% ylabel('deviation (%)');
% hold off

fprintf('Steady-state max displacement for P = 10N: %f\n', yPosLoad);
fprintf('Analytical Euler-Bernoulli max displacement for P = 10N: %f\n', EBdef);
fprintf('Percent Deviation: %f\n', deviation);

