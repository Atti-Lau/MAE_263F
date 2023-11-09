function [q, u, a1, a2] = objfun(q0, u, a1, a2, freeIndex, ...
    tol, refTwist)
global Fg M dt
% Guess
q = q0;
iter = 1; % Number of iterations
error = 10 * tol;
while error > tol
    % Compute reference frame
    [a1Iterate, a2Iterate] = computeTimeParallel(a1, q0, q);
    % Compute reference twist
    tangent = computeTangent(q);
    refTwist_iterate = computeRefTwist(a1, tangent, refTwist);
    % Compute material frame
    theta = q(4:4:end);
    [m1, m2] = computeMaterialDirectors(a1Iterate, a2Iterate, theta);
    % Usual force and Jacobian calculation and Newton's update
    [Fb, Jb] = getFb(q, m1, m2);
    [Ft, Jt] = getFt(q, refTwist_iterate);
    [Fs, Js] = getFs(q);
    Forces = Fb + Ft + Fs + Fg;
    JForces = Jb + Jt + Js;
    % Equations of motion
    f = M / dt * ( (q-q0) / dt - u) - Forces;
    % Jacobian
    J = M / dt^2 - JForces;
    f_free = f(freeIndex);
    J_free = J(freeIndex, freeIndex);
    % Newton's update
    dq_free = J_free \ f_free;
    q(freeIndex) = q(freeIndex) - dq_free;
    % Error
    error = sum( abs( f_free ) );
    fprintf('Iter=%d, error=%f\n', iter, error);
    iter = iter + 1;
end

u = (q - q0) / dt;
a1 = a1Iterate;
a2 = a2Iterate;
end