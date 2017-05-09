function Generate_Mesh_egg(stretch)
% stretch = square of radius

t1=clock();

addpath 'MAT_functions'
addpath 'distmesh'
fd=@(p) p(:,1).^2/stretch+p(:,2).^2/stretch+p(:,3).^2/1-1;
[p,t]=distmeshsurface(fd,@huniform,0.04,1.1*[-1,-1,-1;1,1,1]);
M = find_mass_matrix_surface(t,p);
K = find_stiffness_matrix_surface(t,p);

csvwrite('./p_egg.csv',p)
csvwrite('./tri_egg.csv',t)
[i,j,z] = find(M);
csvwrite('./Egg_M.csv',[i,j,z])
[i,j,z] = find(K);
csvwrite('./Egg_K.csv',[i j z])

t2=clock();
disp(etime(t2,t1));
end