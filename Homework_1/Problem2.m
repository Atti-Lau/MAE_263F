%% Clear
	clear all; close all; clc
    
%% Problem Givens
N = 21; %nodes
dt = .01;
totalTime = 50; %in seconds

%rod params
rodLength = .1;
deltaL = rodLength/(N - 1);
r0 = .001; %rod radius
Young = 1e9;

%sphere node radii
R = zeros(N, 1);
R(:) = deltaL/10;
midNode = (N+1)/2;
R(midNode) = .025;

%densities and other constants
rhoMetal = 7000;
rhoFluid = 1000;
delRho = rhoMetal - rhoFluid;
g = 9.81;
visc = 1000;

%utility vars
ne = N - 1; %num of edges
EI = (Young*pi*r0^4) / 4;
EA = Young*pi*r0^2;

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
    M(2*k-1, 2*k-1) = 4/3*pi*(R(k)^3)*rhoMetal; %assigns on every other entry of diagonal
    M(2*k, 2*k) = M(2*k-1,2*k-1); %mass is repeated in M
end

%viscous damping matrix, C
C = zeros(2*N, 2*N);
for k = 1:N
    C(2*k-1, 2*k-1) = 6*pi*visc*R(k);
    C(2*k, 2*k) = C(2*k-1, 2*k-1);
end

%weight force vector, W
W = zeros(2*N,1);
for k = 1:N
    W(2*k-1) = 0;
    W(2*k) = (-4/3)*pi*(R(k)^3)*delRho*g; %force only in y direction, down is -y
end

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
yPosMid = zeros(nSteps, 1); %stores y-position of mid node
yVelMid = zeros(nSteps, 1); %stores y-velocity of mid node

yPosMid(1) = q(2*midNode); %store initial position
yVelMid(1) = u(2*midNode); %store initial velocity

%% Time Loop

for c = 1:nSteps
    
    fprintf('Time = %f\n', (c-1)*dt);
    
    q = q0; %initial guess
    
    %---Newton-Raphson loop---
    err = 10*tol;
    
    while (err > tol)
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
        
        %viscous terms
        f = f + C*(q - q0)/dt;
        J = J + C/dt;
        
        %weight term
        f = f - W;
        
        %Newton-Raphson calculation
        q = q - J\f;
        
        err = sum(abs(f));
    end %----end of Newton-Raphson loop----
    
    %update old values to new
    u = (q - q0)/dt;
    q0 = q;
    
    %frame-by-frame plotting
%     figure(1);
%     plot( q(1:2:end), q(2:2:end), 'ro-');
%     axis equal
%     drawnow
    
    %store data
    yPosMid(c) = q(4);
    yVelMid(c) = u(4);
    
end

%% Plotting
timeArray = (1:nSteps)*dt;

figure(1);
hold on
title('Final Deformed Shape of Beam')
plot( q(1:2:end), q(2:2:end), 'ro-');
xlabel('x (m)');
ylabel('y (m)');
axis equal
hold off

figure(2)
hold on
title('Position of Middle Node vs. Time (Implicit)')
plot(timeArray, yPosMid)
xlabel('time (s)');
ylabel('y-position (m)');
hold off

figure(3)
hold on
title('Velocity of Middle Node vs. Time (Implicit)')
plot(timeArray, yVelMid)
xlabel('time (s)');
ylabel('y-velocity (m/s)');
hold off

fprintf('Terminal velocity: %f\n', yVelMid(nSteps))
%% Problem 3 Plotting
Narr = [5, 11, 21, 31];
dtArr = [.01, 1, 10, 25];

termVel_Narr = [-.005863, -0.005843, -.005840, -.005840];
termVel_dtArr = [-0.005840, -.005840, -0.005840, -0.005826];

figure(4)
hold on
title('Terminal Velocity vs. Number of Nodes')
plot(Narr, termVel_Narr, 'o-')
xlabel('N');
ylabel('Terminal Velocity (m/s)');
hold off

figure(5)
hold on
title('Terminal Velocity vs. Timestep')
plot(dtArr, termVel_dtArr, 'o-')
xlabel('dt (s)');
ylabel('Terminal Velocity (m/s)');
hold off