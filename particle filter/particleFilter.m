clc
clear all
close all
load a3D;
% we are going to be linear motion to height
% x verically down
% y horizontally right
xMap = x;
yMap = y;
zMap = z;
xRes = xMap(2) - xMap(1);
yRes = yMap(2) - yMap(1);

tend = 150;
dt = 0.1;
time = 0;
N = floor(tend/dt);

px = xMap(floor(0.3*length(xMap)))
py = yMap(floor(0.2*length(yMap)))
pz = max(max(zMap)) + 250

pxs = px;
pys = py;
pzs = pz;

vx =  4;
vy =  8;
vz =  0;

% particle start here
%hnop = 10;
%vnop = 10;

%hmin = -400;
%hmax =  4000;
%vmin = -10;
%vmax =  1000;

%dh = (hmax - hmin)/hnop;
%dv = (vmax - vmin)/vnop;
%k = 1;
%hc = hmin;
%for i=1:hnop
%	hc = hc + dh;
%	vc = vmin; 
%	for j=1:vnop
%		vc = vc + dv;
%		S(k,:) = [hc vc];
%		W(k)   = 1/(hnop*vnop);
%		k = k+1; 
%	end
%end 
%nop = hnop*vnop;

for i=2:N

	ax = 0;
	ay = 0;
	az = 0;
    % dynamics
	px = px + vx*dt;
	vx = vx + ax*dt;
	
	py = py + vy*dt;
	vy = vy + ay*dt;

	pz = pz + vz*dt;
	vz = vz + az*dt;
    
	% sensor
	% -----------------------------------------------

	ix = (px - xMap(1))/xRes;
	iy = (py - yMap(1))/yRes;

	surfaceHeight = 0;
	surfaceHeight = surfaceHeight + (ceil(ix) - ix)*(ceil(iy) - iy)*zMap(floor(ix),floor(iy));
	surfaceHeight = surfaceHeight + (ix - floor(ix))*(ceil(iy) - iy)*zMap(ceil(ix),floor(iy));
	surfaceHeight = surfaceHeight + (ceil(ix) - ix)*(iy - floor(iy))*zMap(floor(ix),ceil(iy));
	surfaceHeight = surfaceHeight + (ix - floor(ix))*(iy - floor(iy))*zMap(ceil(ix),ceil(iy));

	htMeas = pz - surfaceHeight + 0.1*randn(1);

	lhts(i) = surfaceHeight;

	% Resampling 
	%Sout = SUS(S,W);
	%eta = 0;
	%for j=1:length(Sout)
	%	X = Sout(j,:);
	%	Xnext = X + [vhMeas vvMeas]*dt + 2.0*randn(1); % here velocity can be removed
	%	%Zinter = searchMap(hMap,vMap,htmap,Xnext);
        %Zinter  = sMap(localHt,localDs,Xnext);  
        %scale = 0.5;
	%	Wn(j) = exp(-1/(2*100^2)*(scale*((rMeas - htMeas) - (Xnext(2) - Zinter))^2 + (1-scale)*(Xnext(2) - rMeas)^2)); % here exponential can be tried		
	%	eta = eta + Wn(j);
	%	S(j,:) = Xnext;
	%end
	%W = Wn/eta;

	

	pxs(i) = px;
	pys(i) = py;
	pzs(i) = pz;
	htMeasS(i) = htMeas;
	time(i) = time(i-1) + dt;
% 	figure(1),hold on,plot(phs,pvs,'r*'),plot(phs,lhts,'g'),plot(S(:,1),S(:,2),'k.')
%	figure(2),hold on,plot3(S(:,1),S(:,2),W,'b*'),plot(phs,pvs,'r*'),plot(localDs,localHt);grid on
    
end

figure,mesh(xMap,yMap,zMap),hold on,plot3(pxs,pys,pzs,'b*')
figure,plot(htMeasS),grid on


