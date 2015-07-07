V = [0 1; 4 3; 2 5; -4 6];
clf;
hold on;
drawBezier(V);
drawControlVertexes(V);

for i=1:3
    V = increaseGrade(V);
    drawControlVertexes(V);
end

hold off;