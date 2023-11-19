N = 100;
size = 1000;

% sigmaFact = 100;
% mu = sigmaFact/2;
sigmaFact = 100*sqrt(1/12);
mu = 1;
A = mu + sigmaFact*randn(N,N,size);
S = zeros(size,N);
Smax = zeros(size,1);
Sothers = zeros(size,N-1);

for i = 1:size
    S(i,:) = svd(A(:,:,i));
    Smax(i,1) = S(i,1);
    Sothers(i,:) = S(i,2:end);
end

histogram(S(:),'Normalization','pdf', 'FaceColor', "#EDB120")

hold on

sigma = 100*sqrt(1/12);
a = 2 * sigma * sqrt(N);
s = linspace(0,a,300);
f = 4/(pi*a) .* sqrt(1 - s.^2 / a^2);
plot(s,f, "b", LineWidth=2);
title('Normal distribution \mu = 1, \sigma^2 = 100^2/12, N = 100, <S_{max}>_{theory} = 933.3');
ylabel('dw/ds');

sMaxTheory = (mu + sigma^2/(mu * N))*N;
plot([sMaxTheory sMaxTheory], [0 2e-3], "g", LineWidth=2)

