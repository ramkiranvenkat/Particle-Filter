function Sout = SUS(S,W)

N = length(S);
c = W(1);

for i=2:N
	c(i) = c(i-1) + W(i);
end
i = 1;
u = 1/N*rand(1);
for j=1:N
	while(u(j) > c(i))
		i = i+1;
	end
	Sout(j,:) = S(i,:);
    	u(j+1) = u(j) + 1/N;
end
