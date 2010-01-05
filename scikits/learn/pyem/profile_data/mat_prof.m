% Last Change: Mon Jun 04 10:00 AM 2007 J

n   = 1e5;
d   = 30;

x   = randn(n, d);
mu  = randn(n, d);

for i=1:10
    y = gden(x, mu);
end;
