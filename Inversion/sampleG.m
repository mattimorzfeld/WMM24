function out = sampleG(n,a)
y = rand(n,1);
C = 2*(a-1)/sqrt(a);
out = (1/sqrt(a)+(C/2)*y).^2;