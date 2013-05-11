function dy = propVariance(t,y,A)
tempMat = [y(1),0,0,0;...
           0,y(6),0,0;...
           0,0,y(11),0;...
           0,0,0,y(16)];

tempdy = A*tempMat + tempMat*A';
dy = zeros(16,1);
count = 1;
for i = 1:4
    for j = 1:4
        dy(count) = tempdy(i,j);
        count = count + 1;
    end
end
return