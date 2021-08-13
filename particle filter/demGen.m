function [xArray yArray z biasS] = demGen(xSt,ySt,xSp,ySp,nx,ny)

xArray = [xSt:(xSp-xSt)/nx:xSp];
yArray = [ySt:(ySp-ySt)/ny:ySp];

xmid1 = xArray(round(nx/2*rand(1)));
ymid1 = yArray(round(ny/2*rand(1)));
xmid2 = xArray(round(nx/2+nx/2*rand(1)));
ymid2 = yArray(round(ny/2+ny/2*rand(1)));

bias = 0;
for i=1:length(xArray)
	
	if (rem(i,20) == 1)
		ibit = 1;
	else
		ibit = 0;
	end
	for j=1:length(yArray)
		cx = 0.01;
		cy = 0.005;
		
		z(i,j) = 200*exp(-((xArray(i)-xmid1)^2 + (yArray(j)-ymid1)^2)/(2*300^2));
		z(i,j) = z(i,j) + 300*exp(-((xArray(i)-xmid2)^2 + (yArray(j)-ymid2)^2)/(2*600^2));
		z(i,j) = z(i,j)*(sin(cx*xArray(i) + cy*yArray(j)) + 0.05*sin(0.1*sqrt(xArray(i)^2 + yArray(j)^2)) ) + 600 + 5*randn(1);
		if (ibit == 1)||(rem(j,20) == 1)
			bias =  bias + 0.5*randn(1);			
		end
		biasS(i,j) = bias;
		z(i,j) = z(i,j) + bias;
	end 
end

[xm1 ym1 zm1] = mountain((xSp-xSt)/2,(ySp-ySt)/2,350,round(nx/3),round(ny/3));
z(floor(nx/3):floor(nx/3)+length(xm1)-1,floor(nx/3):floor(nx/3)+length(ym1)-1) = z(floor(nx/3):floor(nx/3)+length(xm1)-1,floor(nx/3):floor(nx/3)+length(ym1)-1) + zm1;



