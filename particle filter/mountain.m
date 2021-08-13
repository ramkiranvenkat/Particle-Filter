function [xArray yArray z] = mountain(x,y,maxHt,nx,ny);
xArray = [-x/2:x/nx:x/2];
yArray = [-y/2:y/ny:y/2];

for i=1:length(xArray)
	for j=1:length(yArray)
		z(i,j) = (1 - (xArray(i)^2 + yArray(j)^2)*100/(x^2 + y^2))*maxHt*exp(-((xArray(i))^2 + (yArray(j))^2)/(2*200^2)) + 10*rand(1);
	end
end

 
