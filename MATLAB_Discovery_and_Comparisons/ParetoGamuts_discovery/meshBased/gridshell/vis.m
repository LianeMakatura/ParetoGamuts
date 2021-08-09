g = Gridshell("", 5,5,[1,2]);

vector_1 = [0,0,1] + g.morningSun;
vector_2 = [0,0,1] +  g.eveningSun;
vector_3 = [0,0,1];

figure;
plot3([0 0 0; vector_1(:,1), vector_2(:,1), vector_3(:,1)], ...
      [0 0 0; vector_1(:,2), vector_2(:,2), vector_3(:,2)], ...
      [1 1 1; vector_1(:,3), vector_2(:,3), vector_3(:,3)]);
axis equal; 

b = buffer9;
p = pts9;

visIdx(300, b, p);
visIdx(3000, b, p); 
visIdx(6000, b, p);

function visIdx(idx, b, p)
    figure; 
    side = sqrt(b.rD)+2;
    z = zeros(side);
    z(2:end-1, 2:end-1) = reshape(p(idx, 1:b.rD)*5, sqrt(b.rD), sqrt(b.rD));
    surf(z); 
    title("point " + string(idx)); 
    fprintf("House Theta: " + string(p(idx, b.rD+1)*360) + " deg\n"); 
    fprintf("Morning Perf: " + string(p(idx, b.rD+b.rA+1)) + "\n"); 
    fprintf("Evening Perf: " + string(p(idx, b.rD+b.rA+2)) + "\n")
end