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

for j = 1:Ndetectors
    disp(j)
    myname = ['sig_',num2str(j)];
    %% Load path lengths of detected photons DetS
    filename = sprintf('%s.bin',myname);
    fid = fopen(filename, 'rb');
    Detsig = [fread(fid, 'double')];
    fclose(fid);
    sig_temp = Detsig(1:2:end)+1i*Detsig(2:2:end);
    sig(:,j) = sig_temp;
end
sig = interp1(k_nolin,sig,k_lin);

%% %%%%%%%
%Adding noise
noise_amp = 0.00;
noise1 = raylrnd(noise_amp/1.253*ones(size(sig))) .* exp(1i.*2.*pi.*rand(size(sig)));
noise2 = raylrnd(noise_amp/1.253*ones(size(sig))) .* exp(1i.*2.*pi.*rand(size(sig)));
ref_amp = 1;

%Make the interference with the reference arm and calculate the intensity
I = abs(sig + 1+noise1).^2 - abs(sig - 1+noise2).^2;
lamda0 = (lambda_stop+lambda_start)/2;
delta_k = pi/sqrt(log(2))*((lambda_stop-lambda_start)/2)/(lamda0)^2;
S_k =1;
I = sqrt(S_k).*I;
%% Processing the OCT signal
k = (k_lin);
%Apply low pass filter and hanning window
ksampling = 2*pi/(k(1)-k(2));
if maxDepth == 0
    rawAline = (I).*hanning(length(k));
else
    rawAline = lowpass(I'.*hann(length(k)),maxDepth,ksampling);
end
% rawAline = I;
%Calculate Aline
M10 = length(rawAline(:,1));
OCT = abs(ifft(rawAline,M10));
OCT = OCT(1:floor(length(OCT(:,1))/2)+1,:);
OCT(2:end-1,:) = OCT(2:end-1,:);

%% Displaying image
z = (0:M10/2)/M10.*ksampling;
x = linspace(-radius,radius,Ndetectors);
z=z/2;
figure
imagesc(x,z,10*log10(OCT/max(max(OCT))))
xlabel('Position [cm]')
ylabel('Depth [cm]')
colormap jet

toc
