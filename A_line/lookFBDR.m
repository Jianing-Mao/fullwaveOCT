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
z_min = 2*pi/(k_lin(1)-k_lin(2))/samplePoints;
MM = zeros(samplePoints,samplePoints);
%%
for j = 1:samplePoints
    disp(j)
    myname = ['infi',num2str(j)];
    %% Load path lengths of detected photons DetS
    filename = sprintf('%s_DetS.bin',myname);
    fid = fopen(filename, 'rb');
    DetS = [fread(fid, 'float')];
    fclose(fid);
    %% Load weight of detected photons DetW
    filename = sprintf('%s_DetW.bin',myname);
    fid = fopen(filename, 'rb');
    DetW = [fread(fid, 'float')];
    fclose(fid);
    %% Load likelihood of detected photons DetL
    filename = sprintf('%s_DetL.bin',myname);
    fid = fopen(filename, 'rb');
    DetL = [fread(fid, 'float')];
    fclose(fid);
    
    %% Load the saved photons
    S = DetS(:)';
    W = DetW(:)';
    L=DetL(:)';
    S=S*n_cor;
    
    for i = 1:samplePoints
        MM(i,j) = sum(L(S >= 0+(i-1)*z_min & S < 0+z_min+(i-1)*z_min).*W(S >= 0+(i-1)*z_min & S < 0+z_min+(i-1)*z_min));
    end
    
end
%% %%%%%%%%%%%%%

fun_all = zeros(samplePoints,samplePoints);
fun_all_ori = zeros(samplePoints,samplePoints);
for i=1:1024
    MM70 = MM(i,:);
    p_70 = polyfit([1:1024],MM70,4);
    fun_70 = polyval(p_70,[1:1024]);
    fun_70(fun_70<0)=0;
    fun_70 = fun_70/(sum(fun_70)+eps)*sum(MM70);
    fun_all_ori(i,:) = fun_70;
end
figure;
imagesc(fun_all_ori(1:200,:));colormap jet