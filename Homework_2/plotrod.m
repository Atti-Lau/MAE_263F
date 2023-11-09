function plotrod(q, a1, a2, m1, m2, ctime)

nv = (length(q)+1)/4;
x1 = q(1:4:end);
x2 = q(2:4:end);
x3 = q(3:4:end);

L = sum(sqrt( (x1(2:end) - x1(1:end-1)).^2 + ...
    (x2(2:end) - x2(1:end-1)).^2 + ...
    (x3(2:end) - x3(1:end-1)).^2));

a1 = 0.1*L * a1;
a2 = 0.1*L * a2;
m1 = 0.1*L * m1;
m2 = 0.1*L * m2;
h1 = figure(1);

% set(h1, 'visible', 'off');
clf()
plot3(x1,x2,x3, 'ko-');
hold on
plot3(x1(1),x2(1),x3(1), 'r^');

for c=1:nv-1
    xa = q(4*c-3:4*c-1);
    xb = q(4*c+1:4*c+3);
    xp = (xa+xb)/2;
    p1 = plot3( [xp(1), xp(1) + a1(c,1)], [xp(2), xp(2) + a1(c,2)], ...
        [xp(3), xp(3) + a1(c,3)], 'b--', 'LineWidth', 2);
    p2 = plot3( [xp(1), xp(1) + a2(c,1)], [xp(2), xp(2) + a2(c,2)], ...
        [xp(3), xp(3) + a2(c,3)], 'c--', 'LineWidth', 2);
    p3 = plot3( [xp(1), xp(1) + m1(c,1)], [xp(2), xp(2) + m1(c,2)], ...
        [xp(3), xp(3) + m1(c,3)], 'r-');
    p4 = plot3( [xp(1), xp(1) + m2(c,1)], [xp(2), xp(2) + m2(c,2)], ...
        [xp(3), xp(3) + m2(c,3)], 'g-');
end
hold off

legend([p1,p2,p3,p4], 'a_1','a_2','m_1','m_2');
title(num2str(ctime, 't=%f'));
axis equal
view(3);
xlabel('x');
ylabel('y');
zlabel('z');
drawnow