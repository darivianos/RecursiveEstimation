function dy = propVariance(t,y,A,Q)
tempMat = [y(1),y(5),y(9),y(13);...
           y(2),y(6),y(10),y(14);...
           y(3),y(7),y(11),y(15);...
           y(4),y(8),y(12),y(16)];

tempdy = A*tempMat + tempMat*A' + Q;
dy = zeros(16,1);
count = 1;
for i = 1:4
    for j = 1:4
        dy(count) = tempdy(j,i);
        count = count + 1;
    end
end
return