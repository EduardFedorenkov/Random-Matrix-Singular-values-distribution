N = 100;
size = 1000;

A = -0.5 + rand(N,N,size);
S1 = zeros(size,N);
Smax1 = zeros(size,1);
Sothers1 = zeros(size,N-1);

B = 0 + sqrt(1/12).*randn(N,N,size);
S2 = zeros(size,N);
Smax2 = zeros(size,1);
Sothers2 = zeros(size,N-1);

for i = 1:size
    S1(i,:) = svd(A(:,:,i));
    Smax1(i,1) = S1(i,1);
    Sothers1(i,:) = S1(i,2:end);
    S2(i,:) = svd(B(:,:,i));
    Smax2(i,1) = S2(i,1);
    Sothers2(i,:) = S2(i,2:end);
end

tiledlayout(1,2);

ax1 = nexttile;

histogram(S1(:),'Normalization','pdf', 'FaceColor', "#EDB120")

hold on

sigma = sqrt(1/12);
a = 2 * sigma * sqrt(N);
s = linspace(0,a,300);
f = 4/(pi*a) .* sqrt(1 - s.^2 / a^2);
plot(s,f, "b", LineWidth=2);
title(ax1,'Uniform distribution \mu = 0, \sigma^2 = 1/12, N = 100');
ylabel(ax1,'dw/ds');

ax2 = nexttile;

histogram(S2(:),'Normalization','pdf', 'FaceColor', "#EDB120")

hold on

sigma = sqrt(1/12);
a = 2 * sigma * sqrt(N);
s = linspace(0,a,300);
f = 4/(pi*a) .* sqrt(1 - s.^2 / a^2);
plot(s,f, "b", LineWidth=2);
title(ax2,'Normal distribution \mu = 0, \sigma^2 = 1/12, N = 100');
ylabel(ax2,'dw/ds');
