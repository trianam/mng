function splitBezier(V, t)
%splitBezier Splits the curve in two
%   V are the control vertexes of the curve;
%   t is the parameter of splitting

n = size(V,1) -1;

V1 = [];
V2 = [];

for i=1:n+1
    Q(i,:)=V(i,:);
end

for k=1:n
    V1 = [V1; Q(1,:)];
    for i=1:n+1-k
        Q(i,:) = (1.0-t)*Q(i,:)+t*Q(i+1,:);
    end
end
V1 = [V1; Q(1,:)];

for i=1:n+1
    V2 = [V2; Q(i,:)];
end

disp(V);
disp(V1);
disp(V2);

clf;
hold on;
drawBezier(V);
drawControlVertexes(V);
drawBezier(V1);
drawControlVertexes(V1);
drawBezier(V2);
drawControlVertexes(V2);
hold off;
end

