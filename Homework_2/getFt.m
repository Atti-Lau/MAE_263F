function [Ft, Jt] = getFt(q, refTwist)
global GJ voronoiLength
nv = (length(q)+1) / 4;
Ft = zeros(size(q));
Jt = zeros(length(q), length(q));
for c = 2:nv-1
    node0 = [q(4*c-7), q(4*c-6), q(4*c-5)];
    node1 = [q(4*c-3), q(4*c-2), q(4*c-1)];
    node2 = [q(4*c+1), q(4*c+2), q(4*c+3)];
    theta_e = q(4*c-4);
    theta_f = q(4*c);
    [dF, dJ] = ...
        gradEt_hessEt(node0, node1, node2, ...
        theta_e, theta_f, refTwist(c), ...
        voronoiLength(c), GJ);
    ind = 4*c-7:4*c+3; % Size 11
    Ft(ind) = Ft(ind) - dF;
    Jt(ind, ind) = Jt(ind, ind) - dJ;
end
end