%% Create the files for the simulation
clear
format compact
clc
home
%% load the desired mus and TT parameters of different particle diameters 
load us_4_1300.mat % scattering coefficient of 4 mu_m diameter particles
load us_6_1300.mat % scattering coefficient of 6 mu_m diameter particles
para4 = load('parameters_4.txt'); % TT SPF parameters of 4 mu_m diameter particles
para6 = load('parameters_6.txt'); % TT SPF parameters of 6 mu_m diameter particles

%% set paras
%%% USER CHOICES %%%%%%%%
SAVEON      = 1;        % 1 = save myname_T.bin, myname_H.mci
samplePoints = 1024; %number of wavelength sampling
myname      = 'infi';% name for files: myname_T.bin, myname_H.mci
Nphotons    = 5000;      	% number of photons
Nx          = 200;    	% # of bins in each dimension of cube
Ny          = 200;    	% # of bins in each dimension of cube
Nz          = 200;    	% # of bins in each dimension of cube
binsize     = 0.001;     	% size of each bin, eg. [cm]

% Set Monte Carlo launch flags (not in use)
% Sets position of source with 0 centered on each Aline
xs          = 0;      	% x of source
ys          = 0;        % y of source
zs          = 0.0001; % z of source must start in simulation

% Set detection parameter
radius      = 0.1;     % Half width of the BScan
Ndetectors  = 1;      % Number of Aline per BScan, 1 for Aline simulation

% det_radius  = nm*1e-7*2/pi/atan(beamw/2/flens);    % Width of the beam at the imaging lens
det_radius = 3/2*binsize;

% Cos of the accepted angle
cos_accept = 0.9961;
% manually set launch trajectory
ux0         = 0.0;      % trajectory projected onto x axis
uy0         = 0.0;      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1

% Bias scattering parameter
p  = 0.5;
dx = binsize;
dy = binsize;
dz = binsize;
%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify Monte Carlo parameters
x  = ([0:Nx]'-Nx/2)*dx;
y  = ([0:Ny]'-Ny/2)*dy;
z  = [0:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);
%%%%%%%%%%%%
%% Create Sample
%%% number of mediums
Nt = 9;
%%% set your structure and parameters for each wavelength
for jj=1:samplePoints
    myname = ['infi',num2str(jj)]
    for i=1:Nt
        if i == 1
            nrv(i)   = 1.33;
        elseif i == 4
            nrv(i) = 1.33; % refraction index (set uniform for every medium)
            muav(i) = 0; % absorption coefficient
            musv(i)= us_4(jj); % scattering coefficient
            gv(i) = 0.9; % anisotropy
            gf(i) = para4(jj,1); % forward anisotropy
            gb(i) = para4(jj,2); % backward anisotropy
            alf(i) = para4(jj,3); % forward enhancing
            alb(i) = para4(jj,4); % backward enhancing
            C(i) = para4(jj,5); % balance factor
        elseif i==5
            nrv(i) = 1.33;
            muav(i) = 0;
            musv(i)= us_6(jj);
            gv(i) = 0.9;
            gf(i) = para6(jj,1);
            gb(i) = para6(jj,2);
            alf(i) = para6(jj,3);
            alb(i) = para6(jj,4);
            C(i) = para6(jj,5);
        else
            nrv(i) = 1.4+rand;
            muav(i) = 0;
            musv(i)= 0;
            gv(i) = 0;
            gf(i) = 0;
            gb(i) = 0;
            alf(i) = 0;
            alb(i) = 0;
            C(i) = 0;
        end
    end
    T = double(zeros(Ny,Nx,Nz));
    
    T = T + 1;      %
    
    zsurf = 0.0000;  % position of surface
    %%% set the structure
    for iz=1:Nz % for every depth z(iz)
        
        % water
        if iz<=round(0.01/dz)
            T(:,:,iz) =1;
        end
        
        if iz>round(0.02/dz)
            T(:,:,iz) = 5;
        end
        
    end % iz
    
    
    
    %%
    if SAVEON
        tic
        % convert T to linear array of integer values, v(i)i = 0;
        v = uint8(reshape(T,Ny*Nx*Nz,1));
        
        %% WRITE FILES
        % Write myname_H.mci file
        %   which contains the Monte Carlo simulation parameters
        %   and specifies the tissue optical properties for each tissue type.
        commandwindow
        disp(sprintf('--------create %s --------',myname))
        filename = sprintf('%s_H.mci',myname);
        fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',Nphotons);
        fprintf(fid,'%0.4f\n',p);
        fprintf(fid,'%0.4f\n',Ndetectors);
        fprintf(fid,'%0.6f\n',det_radius);
        fprintf(fid,'%0.6f\n',cos_accept);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
        % launch parameters
        fprintf(fid,'%0.4f\n',xs);
        fprintf(fid,'%0.4f\n',ys);
        fprintf(fid,'%0.4f\n',zs);
        fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
        fprintf(fid,'%0.4f\n',uy0);
        fprintf(fid,'%0.4f\n',uz0);
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',zsurf);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.6f\n',muav(i));
            fprintf(fid,'%0.6f\n',musv(i));
            fprintf(fid,'%0.6f\n',gv(i));
            fprintf(fid,'%0.6f\n',nrv(i));
            fprintf(fid,'%0.6f\n',gf(i));
            fprintf(fid,'%0.6f\n',gb(i));
            fprintf(fid,'%0.6f\n',alf(i));
            fprintf(fid,'%0.6f\n',alb(i));
            fprintf(fid,'%0.6f\n',C(i));
        end
        fclose(fid);
        
        %% write myname_T.bin file
        if jj==1
            filename = sprintf('%s_T.bin',myname);
            disp(['create ' filename])
            fid = fopen(filename,'wb');
            fwrite(fid,v,'uint8');
            fclose(fid);
        end
        toc
    end % SAVEON
end
%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,Ny/2)'; % Tzx

%%
figure('color', 'white'); clf
sz = 12;  fz = 10;
imagesc(x,z,Tzx,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('\rm Tissue')
colorbar
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%%


text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([xmin xmax zmin zmax])

savefig(strcat(myname, '.fig'));
disp('done')


