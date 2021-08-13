function z = searchMap(x,y,ht,Bel)
xlen = length(x);
xinit = x(1);
i = 2;
while (x(i) == xinit)
	i = i + 1;
end
ylen = i-1;

xAlone = x(1:ylen:xlen);
yAlone = y(1:ylen);

xlen = length(xAlone);
ylen = length(yAlone);

dx = (xAlone(xlen) - xAlone(1))/xlen;
dy = (yAlone(ylen) - yAlone(1))/ylen;

iXSearch = (Bel(1)-xAlone(1))/dx;
iYSearch = (Bel(2)-yAlone(1))/dy;

if ((iXSearch < 0) || (iYSearch < 0))
display('Out of Input Map Range');
z = 0;
else
itr = round(iXSearch)*ylen + round(iYSearch);	
z = ht(itr);
end

