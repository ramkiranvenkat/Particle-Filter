function z = sMap(ht,ds,bel)
DsRes   =  0.1;
i = (bel(1) - ds(1))/DsRes;
height = (ceil(i) - i)*ht(floor(i)) + (i-floor(i))*ht(ceil(i));
z = bel(2) - height;