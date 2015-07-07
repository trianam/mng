function drawBezier(V)
%drawBezier draw the Bezier curve
%   V is the vector of control vertex
%       (matrix n X 2: column 1 = x; column 2 = y)
n = size(V,1) -1;

C = [];
for t = 0:0.01:1
    C = [C; deCasteljau(V, n, t)];
end

if (size(C,2) == 2)
    plot(C(:,1), C(:,2));
elseif (size(C,2) == 3)
    plot3(C(:,1), C(:,2), C(:,3));
else
    disp('ERR');
end

end

