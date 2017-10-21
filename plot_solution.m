clear all;

load 'cells.m';
load 'solution005.m';
u=solution005;

figure(1)
trimesh(cells,u(:,1),u(:,2),u(:,3));
axis square;
axis tight;
xlabel('x');
ylabel('y');
zlabel('\rho');


figure(2)
trimesh(cells,u(:,1),u(:,2),u(:,4));
axis square;
axis tight;
xlabel('x');
ylabel('y');
zlabel('\rho u');

figure(3)
trimesh(cells,u(:,1),u(:,2),u(:,5));
axis square;
axis tight;
xlabel('x');
ylabel('y');
zlabel('\rho v');

figure(4)
trimesh(cells,u(:,1),u(:,2),u(:,6));
axis square;
axis tight;
xlabel('x');
ylabel('y');
zlabel('E');