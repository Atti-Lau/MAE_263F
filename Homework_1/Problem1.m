%% Clear
	clear all %#ok<*CLALL>
	close all
    clc

%% Node and DOF, Problem Setup
runExplicit = true; %set this flag to run explicit solver

%---Algorithm params---
N = 3; %# of nodes
nDOF = N*2; %for 2D

dtIm = .1; %implicit timestep
dtEx = 1e-4; %explicit timestep

totalTime = 10; %in seconds
%---problem parameters---
%rod params
rodLength = .1;
deltaL = rodLength/(N-1); %space between spheres
rodR = .001;
YoungsMod = 1e9;
grav = 9.81;

%sphere radii
r1 = .005;
r2 = .025;
r3 = .005;

%material densities
rhoMetal = 7000; 
rhoFluid = 1000;

delRho = rhoMetal - rhoFluid;

visc = 1000;

%Utility parameters
numEdges = N - 1; %number of edges
EI = YoungsMod*pi*rodR^4 / 4; %rod bending stiffness
EA = YoungsMod*pi*rodR^2; %Newton WHAT?

%% Algorithm Setup
%-----Geometry, initial positioning----
nodes = zeros(N, 2); %each row is a sphere's x and y

for c = 1:N %sets up initial positions of nodes
    nodes(c,1) = (c-1)*deltaL; %x-positions of nodes
    nodes(c,2) = 0; %all start at top, so y = 0
end

%mass matrix M
M = zeros(nDOF, nDOF);

M(1,1) = (4/3)*pi*(r1^3)*rhoMetal; %sphere 1
M(2,2) = M(1,1);
M(3,3) = (4/3)*pi*(r2^3)*rhoMetal; %sphere 2
M(4,4) = M(3,3);
M(5,5) = (4/3)*pi*(r3^3)*rhoMetal; %sphere 3
M(6,6) = M(5,5);

%damping matrix C
C = zeros(nDOF, nDOF);

C(1,1) = 6*pi*visc*r1; %sphere 1
C(2,2) = C(1,1);
C(3,3) = 6*pi*visc*r2; %sphere 2
C(4,4) = C(3,3);
C(5,5) = 6*pi*visc*r3; %sphere 3
C(6,6) = C(5,5);

%weight vector, W
W = zeros(nDOF, 1);
W(2) = -(4/3)*pi*(r1^3)*(delRho)*grav;
W(4) = -(4/3)*pi*(r2^3)*(delRho)*grav;
W(6) = -(4/3)*pi*(r3^3)*(delRho)*grav;

%Initialize DOFs
q0 = zeros(nDOF, 1);

for c = 1:N %assign old/initial DOF vector
    q0(2*c - 1) = nodes(c,1);
    q0(2*c) = nodes(c,2);
end

u0 = zeros(nDOF,1); %initial/old velocity

tol = EI/(rodLength^2) * 1e-3; %small tolerance

%% Time Algorithm Loop (Newton-Raphson, Implicit)

nSteps = round(totalTime/dtIm + 1);

figure(1); %for sphere display
hold on;

%store y-position and velocity of R2:
r2pos = zeros(nSteps, 1);
r2vel = zeros(nSteps, 1);

for c = 1:nSteps %loop to find position at each timestep
    
    currentTime = (c-1)*dtIm;
    fprintf('Time = %f\n', currentTime);
    
    %guess
    q = q0; %set to old DOF
    
    %Newton-Raphson loop
    err = 10*tol;
    while err > tol
        %--total force term, M*q double-dot--
        f = (M/dtIm) * ((q - q0)/dtIm - u0); %first term, net force
        J = M / (dtIm^2); %Jacobian for f, the inertia one in the slides
        
        %elastic forces
        
        %---bending spring forces, grad elastic wrt q---
        %linear between nodes 1 and 2
        xk = q(1);
        yk = q(2);
        xkp1 = q(3);
        ykp1 = q(4);
        l_k = deltaL;
        dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
        
        f(1:4) = f(1:4) + dF;
        J(1:4, 1:4) = J(1:4, 1:4) + dJ;
        
        %linear between nodes 2 and 3
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        l_k = deltaL;
        dF = gradEs(xk, yk, xkp1, ykp1, l_k, EA);
        dJ = hessEs(xk, yk, xkp1, ykp1, l_k, EA);
        
        f(3:6) = f(3:6) + dF;
        J(3:6, 3:6) = J(3:6, 3:6) + dJ;
        
        %bending spring force at node 2
        xkm1 = q(1);
        ykm1 = q(2);
        xk = q(3);
        yk = q(4);
        xkp1 = q(5);
        ykp1 = q(6);
        curvature0 = 0; %no initial curvature
        l_k = deltaL;
        dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
        dJ = hessEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, l_k, EI);
        
        f(1:6) = f(1:6) + dF;
        J(1:6, 1:6) = J(1:6, 1:6) + dJ;
        
        %--viscous forces--
        f = f + C*(q-q0)/dtIm;
        J = J + C/dtIm;
        
        %weight
        f = f - W;
        
        %Update old values
        q = q - J \ f; %the big Newton-Raphson calculation
        err = sum(abs(f));
        
    end
    
    %calc new velocity
    u = (q - q0)/dtIm;
    
    %store info for plotting
    r2pos(c) = q(4);
    r2vel(c) = u(4);
    
    %----plotting-----
    %frame by frame
%     figure(1);
%     plot(q(1:2:end), q(2:2:end), 'ro-');
%     axis equal
%     xlabel('x (m)');
%     ylabel('y (m)');
%     drawnow
    
    %assignment 1, time steps, implicit
    if (ismember(currentTime, [0 .01 .05 .1 1 10]))
        figure(1)
        labelName = ['t = ', num2str(currentTime), 's'];
        plot(q(1:2:end), q(2:2:end), 'o-', 'DisplayName', labelName);
        axis equal
        xlabel('x (m)');
        ylabel('y (m)');
    end
    
    %----update old params---
    q0 = q;
    u0 = u;
end

title('Position of Nodes at Select Times (Implicit)')
legend();


%% Implicit Plotting
timeArrayIm = (1:nSteps)*dtIm;

%Plotting for R2

figure(2)
title('Position of R2 vs. Time (Implicit)')
plot(timeArrayIm, r2pos)
xlabel('time (s)');
ylabel('y-position of R2 (m)');

figure(3)
title('Velocity of R2 vs. Time (Implicit)')
plot(timeArrayIm, r2vel)
xlabel('time (s)');
ylabel('y-velocity of R2 (m/s)');

%% Time Loop (Explicit)
if (runExplicit)
    %re-initialize
    q0 = zeros(nDOF, 1);
    u0 = zeros(nDOF, 1);
    
    for c = 1:N %assign old/initial DOF vector
        q0(2*c - 1) = nodes(c,1);
        q0(2*c) = nodes(c,2);
    end
    
    q = q0;
    u = (q - q0)/dtEx; %initial/old velocity
    
    nSteps = round(totalTime/dtEx + 1);
    
    figure(4); %for sphere display
    hold on;
    
    %store y-position and velocity of R2:
    r2pos = zeros(nSteps, 1);
    r2vel = zeros(nSteps, 1);
    
    for c = 1:nSteps %loop to find position at each timestep CHANGE TO NUMSTEPS IN FINAL
        
        currentTime = (c-1)*dtEx;
        fprintf('Time = %f\n', currentTime);
        
        %q = q0;
        %u = u0;
        
        potE = zeros(nDOF, 1);
        
        %---bending spring forces, grad elastic wrt q---
        %linear between nodes 1 and 2
        for k = 1:(N-1)
            xk = q(2*k - 1);
            yk = q(2*k);
            xkp1 = q(2*k + 1);
            ykp1 = q(2*k + 2);
            dF = gradEs(xk, yk, xkp1, ykp1, deltaL, EA);
            ind = [2*k - 1; 2*k; 2*k + 1; 2*k + 2]; %sets to correct 4x4 spot in matrix
            potE(ind) = potE(ind) + dF;
        end
        
        %bending spring force at node 2
        for k = 2:(N-1)
            xkm1 = q(2*k - 3);
            ykm1 = q(2*k - 2);
            xk = q(2*k - 1);
            yk = q(2*k);
            xkp1 = q(2*k + 1);
            ykp1 = q(2*k + 2);
            curvature0 = 0;
            dF = gradEb(xkm1, ykm1, xk, yk, xkp1, ykp1, curvature0, deltaL, EI);
            ind = [2*k - 3; 2*k - 2; 2*k - 1; 2*k; 2*k + 1; 2*k + 2];
            potE(ind) = potE(ind) + dF;
        end
        
        %-----Explicit Equation----
        overM = zeros(nDOF, nDOF);
        for k = 1:nDOF
            overM(k,k) = 1/(M(k,k));
        end
        
        q = q0 + dtEx*(u0 - (dtEx*overM)*(potE + C*u0 - W)); %Something wrong here
        
        %calc new params, update
        u0 = (q - q0)/dtEx;
        q0 = q;
        
        %store info for plotting
        r2pos(c) = q(4);
        r2vel(c) = u0(4);
        
        %----plotting-----
        %frame by frame
%         figure(4);
%         plot(q(1:2:end), q(2:2:end), 'ro-');
%         axis equal
%         xlabel('x (m)');
%         ylabel('y (m)');
%         axis equal
%         drawnow
        
        %assignment 1, time steps, explicit
        if (ismember(currentTime, [0 .01 .05 .1 1 10]))
            figure(4)
            labelName = ['t = ', num2str(currentTime), 's'];
            plot(q(1:2:end), q(2:2:end), 'o-', 'DisplayName', labelName);
            axis equal
            xlabel('x (m)');
            ylabel('y (m)');
        end
    end
    
    legend();

end %end of explicit solver

%Plotting
timeArrayEx = (1:nSteps)*dtEx;

figure(5)
title('Position of R2 vs. Time (Explicit)')
plot(timeArrayEx, r2pos)
xlabel('time (s)');
ylabel('y-position of R2 (m)');

figure(6)
title('Velocity of R2 vs. Time (Explicit)')
plot(timeArrayEx, r2vel)
xlabel('time (s)');
ylabel('y-velocity of R2 (m/s)');
