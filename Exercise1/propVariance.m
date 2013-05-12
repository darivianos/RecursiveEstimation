function dy = propVariance(t,y,A,Q,L)
tempMat = [y(1),y(2),y(3),y(4);...
           y(5),y(6),y(7),y(8);...
           y(9),y(10),y(11),y(12);...
           y(13),y(14),y(15),y(16)];

tempdy = A*tempMat + tempMat*A' + L*Q*L';
dy = zeros(16,1);
count = 1;
for i = 1:4
    for j = 1:4
        dy(count) = tempdy(i,j);
        count = count + 1;
    end
end
return