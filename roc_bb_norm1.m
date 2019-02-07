%%%%%%%INPUT: 

%%%%%%%nx: number of observations from the first class
%%%%%%%ny: number of observations from the second class
%%%%%%%nz: number of observations from the third class

%%%%%%%x: observations from the first class
%%%%%%%y: observations from the second class
%%%%%%%z: observations from the third class
%%%%%%%rep: resample size
%%%%%%%grid: the length of grid points

% Example:
nx=50;
ny=50;
nz=50;

%(-4, 3, 3, 4, 5, 5)
x=normrnd(-4,3,nx,1);
y=normrnd(3, 4, ny, 1);
z=normrnd(5, 5, nz, 1);

grid=0.0005
rep=100000

t= [grid:grid:1-grid] % FPF vector
ot=ones(length(t),1); % vector of 1 with the same length as vector t
onx=ones(nx,1);ony=ones(ny,1);onz=ones(nz,1);

%%%%%%%%%%%%%%%%%%%BB estimate of ROC, AUC;
for run=1:rep
    % note: to generate Dirichlet weight vectors p and q
p=exprnd(1,1,ny);p=p/sum(p);
q=exprnd(1,1,nx);q=q/sum(q);
r=exprnd(1,1,nz);r=r/sum(r);

s=p*(y*onx'< ony*x');
r2(run,:) = q*(s'*ot'<onx*t);

v=p*(y*onz'< ony*z');
r3(run,:) = r*(v'*ot'<onz*t);

in(run,:)=r2(run,:).*(1-r3(run,:));

vus(run)=auc(in(run,:), grid)

end;

r2bb = mean(r2);
r3bb = mean(r3);
inbb=r2bb.*(1-r3bb);

vusbb=auc(inbb, grid);

vusbb1=sum(inbb)*grid;

