image = phantom(64);
x = linspace(-5,5,64);
[xgrid,ygrid] = ndgrid(x,x);
image = image.* exp(1i * (xgrid.^2 + ygrid.^2));

tic;
unwrapped = phaseUnwrap(angle(image));
toc;

figure(11);
imagesc(unwrapped);
axis square;