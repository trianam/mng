function Vs = increaseGrade(V)
%increaseGrade increase the grade of the control poligon
%   V is the vector of control vertex
%       (matrix n X 2: column 1 = x; column 2 = y)

n = size(V,1);

%{
W = [];
W = [W; 0 0];
for i=1:n
    W = [W; V(i,:)];
end
W = [W; 0 0];
%}

W = [0 0; V; 0 0];
Vs = [];
for i=2:n+2
    Vs = [Vs; ((i-2)/(n)*W(i-1,:) + (1-((i-2)/n))*W(i,:))]; 
end

end
