function out = gden(x, mu)

% Last Change: Mon Jun 04 10:00 AM 2007 J
[n, d] = size(x);
[nm, dm] = size(mu);
if nm ~= n
    out = sum(x-repmat(mu, n, 1), 1);
else
    out = sum(x-mu, 1);
end;
