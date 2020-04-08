n = 64;
image = ones(n,n,n);%phantom3d(n);
x = linspace(-1,1,n);
[xgrid,ygrid,zgrid] = ndgrid(x,x,x);
image = image.* exp(1i * (xgrid.^2 + ygrid.^2 + zgrid.^2));
image(floor(n/4):floor(3*n/4),floor(n/4):floor(3*n/4),floor(n/4):floor(3*n/4)) = 0;
%%
tic;
unwrapped = phaseUnwrap(angle(image));
toc;

figure(11);
imagesc(unwrapped(:,:,32));
axis square;