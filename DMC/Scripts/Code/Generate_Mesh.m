addpath ./MAT_functions/
addpath ./distmesh/

fd=@(p) dsphere(p,0,0,0,1);
[p,t]=distmeshsurface(fd,@huniform,0.08,1.1*[-1,-1,-1;1,1,1]);
M = find_mass_matrix_surface(t,p);
K = find_stiffness_matrix_surface(t,p);

csvwrite('./p.csv',p)
csvwrite('./tri.csv',t)
[i,j,z] = find(M);
csvwrite('./Sphere1_M.csv',[i,j,z])
[i,j,z] = find(K);
csvwrite('./Sphere1_K.csv',[i j z])

figure
x = csvread('./xout.csv');
h = trimesh(t,p(:,1),p(:,2),p(:,3));
set(h,'EdgeColor','k');
hold on
trisurf(t,p(:,1),p(:,2),p(:,3),x);
set(h,'EdgeColor','k');
colormap pink
shading interp
caxis([-2,2]); view(3)

% y = csvread('./yout.csv');
% plot3(y(:,1),y(:,2),y(:,3),'go')

figure
xmean = csvread('./xmean.csv');
h = trimesh(t,p(:,1),p(:,2),p(:,3));
set(h,'EdgeColor','k');
hold on
trisurf(t,p(:,1),p(:,2),p(:,3),xmean);
set(h,'EdgeColor','k');
colormap bone
shading interp
caxis([-2,2]); view(3)
