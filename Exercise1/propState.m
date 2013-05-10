function dy = propState(t,y,uv,ur,B)
dy = zeros(3,1);    % a column vector
dy(1) = y(4) * uv * cos(ur) * cos(y(3));
dy(2) = y(4) * uv * cos(ur) * sin(y(3));
dy(3) = -(y(4) * uv * sin(ur))/B;
dy(4) = 0;
return