function [Fb, Jb] = getFb(q, m1, m2)
global kappaBar EI voronoiLength
nv = (length(q)+1) / 4;
Fb = zeros(size(q));
Jb = zeros(length(q), length(q));
for c = 2:nv-1
    node0 = [q(4*c-7), q(4*c-6), q(4*c-5)];
    node1 = [q(4*c-3), q(4*c-2), q(4*c-1)];
    node2 = [q(4*c+1), q(4*c+2), q(4*c+3)];
    m1e = m1(c-1,:);
    m2e = m2(c-1,:);
    m1f = m1(c,:);
    m2f = m2(c,:);
    [dF, dJ] = ...
        gradEb_hessEb(node0, node1, node2, ...
        m1e, m2e, m1f, m2f, ...
        kappaBar(c,:), voronoiLength(c), EI);
    ind = 4*c-7:4*c+3; % Size 11
    Fb(ind) = Fb(ind) - dF;
    Jb(ind, ind) = Jb(ind, ind) - dJ;
end
end