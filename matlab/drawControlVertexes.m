function drawControlVertexes(V)
%drawControlVertexes draw the control vertexes
%   V is the vector of control vertex
%       (matrix n X 2: column 1 = x; column 2 = y)

%{

n = size(V,1) -1;

Vx=[];
Vy=[];

for i=1:n+1
    Vx = [Vx V(i,1)];
    Vy = [Vy V(i,2)];
end

plot(Vx, Vy);
%}

if (size(V,2) == 2)
    plot(V(:,1), V(:,2));
elseif (size(V,2) == 3)
    plot3(V(:,1), V(:,2), V(:,3));
else
    disp('ERR');
end

end