% clear all
tic
myname = 'infi1';
eps = 1e-10;

filename = sprintf('%s_H.mci',myname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% unit: cm
lambda_start =1250e-7;
lambda_stop = 1350e-7;
samplePoints=1024;
n_cor = 1.33;
lambda_lin = linspace(lambda_start,lambda_stop,samplePoints);
k_nolin = 2*pi./lambda_lin;
k_lin = linspace(2*pi/lambda_start,2*pi/lambda_stop,samplePoints);
sig = zeros(samplePoints,1);
eps=1e-10;
%% parameters
n = 1;
Nphotons = A(n); n = n + 1;
p = A(n); n = n + 1;
Ndetectors = A(n); n = n + 1;
det_radius = A(n); n = n + 1;
cos_accept = A(n); n = n + 1;
Nx = A(n); n = n + 1;
Ny = A(n); n = n + 1;
Nz = A(n); n = n + 1;
dx = A(n); n = n + 1;
dy = A(n); n = n + 1;
dz = A(n); n = n + 1;
xs = A(n); n = n + 1;
ys = A(n); n = n + 1;
zs = A(n); n = n + 1;
ux0 = A(n); n = n + 1;
uy0 = A(n); n = n + 1;
uz0 = A(n); n = n + 1;
radius = A(n); n = n + 1;
zsurf = A(n); n = n + 1;
Nt = A(n);  n = n + 1;

BDR_all = zeros(samplePoints,samplePoints,Ndetectors);
BDR = zeros(samplePoints,samplePoints);
for j = 1:Ndetectors
    disp(j)
    jj=j;
    myname = ['BDR_',num2str(jj)];
    %% Load path lengths of detected photons DetS
    filename = sprintf('%s.bin',myname);
    fid = fopen(filename, 'rb');
    DetBDR = [fread(fid, 'double')];
    fclose(fid);
    
    DetBDR = reshape(DetBDR,samplePoints,[])';
    BDR_all(:,:,j) = DetBDR;
    BDR(:,j) =DetBDR(:,Ndetectors);
end
figure;imagesc(10*log10(BDR));colormap jet