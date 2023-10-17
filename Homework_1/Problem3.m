%% Clear
	clear all; close all; clc
    
%% Problem Givens
N = 50; %nodes
dt = .01;
totalTime = 1; %in seconds

%rod params
rodLength = 1;
deltaL = rodLength/(N - 1);
R = .013;
r = .011; %rod radius
Young = 70e9;
rho = 2700;

PLoad = 20000; %point load, in N
d = .75; %location of point load, +x

%utility vars
ne = N - 1; %num of edges
EI = Young*(pi/4)*(R^4 - r^4);
EA = Young*pi*(R^2 - r^2);

%% Vector Setup
%node matrix, each row is initial x and y of a node
nodes = zeros(N, 2);
for k = 1:N
    nodes(k,1) = (k-1)*deltaL; %initial x of node
    nodes(k,2) = 0; %initial y of node, all = 0
end

%mass matrix M
M = zeros(2*N, 2*N);
for k = 1:N
    M(2*k-1, 2*k-1) = (pi*(R^2 - r^2)*rodLength*rho)/(N-1); %assigns on every other entry of diagonal
    M(2*k, 2*k) = M(2*k-1,2*k-1); %mass is repeated in M
end

%point force vector, P
P = zeros(2*N,1);
PNode = round(d/deltaL);
P(2*PNode) = -PLoad;

%initial DOF vector
q0 = zeros(2*N, 1);
for k = 1:N
    q0(2*k - 1) = nodes(k,1); %x-coord
    q0(2*k) = nodes(k,2); %y-coord
end

%set new position and velocity
q = q0;
u = (q - q0)/dt;

tol = (EI/(rodLength^2))*1e-3;
nSteps = round(totalTime/dt);

%vectors to store data over time
yPosMax = zeros(nSteps, 1); %stores max y-position(columns) over time (rows)
yPos(1) = q(2); %store initial position

%% Time Loop

for c = 1:nSteps
    
    fprintf('Time = %f\n', (c-1)*dt);
    
    %boundary conditions
    fixedDOF = [1; 2; 2*N];
    freeDOF = 3:(2*N - 1);
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
    %q(freeInd) = qFree;
    %q(fixInd) = [0, 0, 0];
    
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

%% Plotting

timeArr = (1:nSteps)*dt;

% maxYs = zeros(nSteps, 1);
% for c = 1:nSteps
%     maxYs(c) = max(yPos(c,:));
% end

figure(2)
hold on
title('Absolute Maximum Vertical Displacement vs. Time')
plot(timeArr, yPosMax)
xlabel('time (s)');
ylabel('maximum y-displacement (m)');
hold off

fprintf('Steady-state max displacement: %f\n', yPosMax(nSteps));
