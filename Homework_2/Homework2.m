%% Clear
clear all; close all; clc
%% Global Vars
global Fg M dt
global kappaBar EI voronoiLength % bending stiffness
global GJ % twisting stiffness
global EA refLen % stretching stiffness

%% Problem Inputs

%node variables
N = 30; % number of nodes
Ne = N - 1; % number of edges
ndof = 3*N + Ne; % number of DOF = the nodes and edges = 4*N - 1
dt = 0.01; % Timestep

% Geometry params
rodLength = 0.2;
natR = 0.02; % natural curvature radius
r0 = 0.001; % cross-sectional radius

% material params
rho = 1000; % Density
Young = 10e6; % Young's Modulus
G = Young/3; % Shear modulus

g = [0; 0; -9.81]; %gravity

totalTime = 5; % total simulation time
nsteps = round(totalTime/dt);

% stiffness variables
EI = (Young*pi*r0^4)/4; % Bending stiffness
GJ = (G*pi*r0^4)/2; % Shearing stiffness
EA = Young*pi*r0^2; % Stretching stiffness

tol = EI/(rodLength^2) * 1e-6; %tolerance

%% Vector/Matrix Setups

%initial geometry vector
nodes = zeros(N, 3);
dTheta = (rodLength/natR)*(1/Ne);
for k = 1:N %initial geometry of rod with curvature
    nodes(k,1) = natR*cos((k-1)*dTheta);
    nodes(k,2) = natR*sin((k-1)*dTheta);
    nodes(k,3) = 0;
end

% Initial DOF vector
q0 = zeros(ndof, 1);
for k=1:N
    ind = [4*k-3; 4*k-2; 4*k-1];
    q0(ind) = nodes(k,:); %as nodes has 3 columns, this maps to above ind
end

% Mass matrix
totalM = pi*(r0^2)*rodLength*rho; %total mass
dm = totalM/Ne; % mass per edge
massVector = zeros(ndof, 1);
for k = 1:N % Loop over nodes
    ind = [4*k-3; 4*k-2; 4*k-1]; %x, y, z indices
    if k == 1
        massVector(ind) = dm/2; %first and last nodes only have half mass
    elseif k == N
        massVector(ind) = dm/2;
    else
        massVector(ind) = dm;
    end
end
for k=1:Ne
    massVector(4*k) = (1/2)*dm*r0^2;
end
M = diag(massVector); % ndof x ndof-sized mass matrix

% reference lengths (edge lengths)
refLen = zeros(Ne, 1);
for k = 1:Ne % loop over the edges
    dx = nodes(k+1,:) - nodes(k,:); %just distances between nodes
    refLen(k) = norm(dx);
end

% Voronoi lengths
voronoiLength = zeros(N, 1); %node-dependent, used for bending and twisting
for k=1:N % loop over the nodes
    if k==1
        voronoiLength(k) = 0.5*refLen(k);
    elseif k==N
        voronoiLength(k) = 0.5*refLen(k-1);
    else
        voronoiLength(k) = 0.5*(refLen(k-1) + refLen(k)); %average of surrounding refLen's
    end
end

%Weight force matrix
Fg = zeros(ndof, 1);
for k=1:N % loop over the nodes
    ind = [4*k-3; 4*k-2; 4*k-1];
    Fg(ind) = massVector(ind).*g;
end

%---Initial reference frame (using space parallel transport)---
a1 = zeros(Ne, 3); %first reference director
a2 = zeros(Ne, 3); %second reference director
tangent = computeTangent(q0); %outputs tangent of all edges in 3 directions

%computing reference directors for first edge
t0 = tangent(1,:); % first tangent edge
t1 = [0;0;-1]; % some vector
a1Tmp = cross(t0, t1); % perpendicular to t0 and t1
if abs(a1Tmp) < 1e-6
    t1 = [0;1;0];
    a1Tmp = cross(t0, t1);
end
a1(1,:) = a1Tmp/norm(a1Tmp);
a2(1,:) = cross(tangent(1,:), a1(1,:));

%reference directors for second edge onward
for k=2:Ne
    t0 = tangent(k-1,:); % tanget on (k-1) edge
    t1 = tangent(k,:); % tanget on k edge
    a1_0 = a1(k-1, :);
    a1_l = parallel_transport(a1_0, t0, t1);
    a1(k,:) = a1_l/norm(a1_l);
    a2(k,:) = cross(t1, a1(k,:));
end

% material frame
theta = q0(4:4:end);
[m1, m2] = computeMaterialDirectors(a1, a2, theta);

% reference twist
refTwist = zeros(N, 1); %kept 0 as a1 and a2 were initialized in parallel transport step

% Natural curvature
kappaBar = getkappa(q0, m1, m2);

% Fixed and free DOFs
fixedIndex = 1:7;
freeIndex = 8:ndof;


%% Timestep Solving

Nsteps = round(totalTime/dt);
ctime = 0; % current time
endZ = zeros(Nsteps, 1); % z-coordinate of last node

% Initialize position and velocity
q = q0;
u = zeros(size(q));

for timeStep = 1:Nsteps
    fprintf('Current time=%f\n', ctime);
    
    [q, u, a1, a2] = objfun(q0, u, a1, a2, freeIndex, ...
        tol, refTwist);
    ctime = ctime + dt;
    
    q0 = q; % Update q
    
    endZ(timeStep) = q(end); % Store z-coordinate at this timestep
    
    if mod(timeStep, 100) == 0
        theta = q(4:4:end);
        [m1, m2] = computeMaterialDirectors(a1, a2, theta);
        plotrod(q, a1, a2, m1, m2, ctime);
    end
end

%% Visualization
figure(2);
timearray = (1:1:Nsteps) * dt;
plot(timearray, endZ, 'ro-');
box on
xlabel('Time, t [sec]');
ylabel('z-coord of last node, \delta_z [m]');


totalM = pi * r0^2 * rodLength * rho;
dm = totalM/Ne; % mass per edge
massVector = zeros(ndof, 1);
for c = 1:N % Loop over nodes
    ind = [4*c-3; 4*c-2; 4*c-1];
end