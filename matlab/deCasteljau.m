function C = deCasteljau(V, n, t)
%deCasteljau calculate the value of the Bezier curve in a point
%   n is the grade of the curve;
%   V is the vector of the n+1 control points
%       (matrix n+1 X 2: column 1 = x; column 2 = y);
%   t is the parameter on wich I want to evaluate the curve
for i=1:n+1
    Q(i,:)=V(i,:);
    %{
    disp('V');
    disp(V(i,:));
    disp('Q');
    disp(Q(i,:));
    %}
end

for k=1:n
    for i=1:n+1-k
        Q(i,:) = (1.0-t)*Q(i,:)+t*Q(i+1,:);
    end
end
C = Q(1,:);
end

