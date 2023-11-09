function refTwist = computeRefTwist(a1, tangent, refTwist)
[ne, ~] = size(a1);

for c = 2 : ne % internal nodes
    u0 = a1(c-1,:); % a1 vector of the previous edge
    u1 = a1(c,:); % a1 of the current edge
    t0 = tangent(c-1,:); % tangent of previous edge
    t1 = tangent(c, :); % tangent of current edge
    ut = parallel_transport(u0, t0, t1);
    % Are ut and u1 the same? NO!
    % Method 1: OK but not so safe?
    % refTwist(c) = signedAngle(ut, u1, t1); % On BruinLearn
    % Method 2
    ut = rotateAxisAngle( ut, t1, refTwist(c));
    refTwist(c) = refTwist(c) + signedAngle(ut, u1, t1);
end

end

function vNew = rotateAxisAngle( v, z, theta )
if (theta == 0)
    vNew = v;
else
    c = cos(theta);
    s = sin(theta);
    vNew = c*v + s*cross(z,v) + dot(z,v)*(1-c)*z;
end

end