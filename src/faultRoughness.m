% CERI 7104/8104 Data Analysis in Geophysics
%
%
% fault roughness analysis
%
% @author thgoebel, CERI U of Memphis
clear all

%=======================1==================================================
%                  parameters
%==========================================================================
data_dir = "~/Documents/teaching/CERI/DataAnalysis_7104/data/";
file_in  = "fault_roughness.txt";

plot_surf = 0;
zMin = -40; % only for surface plots
zMax =  40;

%=======================2==================================================
%                   load data, gridding, detrending etc.
%==========================================================================
XYZ = load(strcat(data_dir,file_in));
%A% remove nans
XYZ =  XYZ(~isnan(XYZ(:,3)),:);
Z     =  XYZ(:,3);
X     =  XYZ(:,1);
Y     =  XYZ(:,2);
a_x_uni = unique( X);
a_y_uni = unique( Y);
sprintf( '%s %.1f %.1f','range of topography:',min(Z), max(Z))
% spatial increments
dx  = X(2) - X(1);
fprintf( '%s%.2f', 'grid spacing in x :', dx)
dy  = a_y_uni(2) - a_y_uni(1);
fprintf( '%s%.2f', 'grid spacing in y :', dy)
%B% interpolation, fill gaps
[XX, YY] = meshgrid( a_x_uni, a_y_uni);
ZZ       = griddata(X, Y, Z, XX, YY, 'linear');

%D% detrend
ZZ_detr = detrend_2d(ZZ);

%C%visualize raw data
if plot_surf == 1
    figure(1)
    subplot(221)
    pcolor( XX, YY, ZZ); 
    title( 'raw data')
    shading interp;
    cbar = colorbar('eastoutside');
    cbar.Label.String = 'Roughness';
    
    subplot(222)
    title( 'detrended')
    pcolor( XX, YY, ZZ_detr); 
    colormap hot
    shading flat;
    cbar = colorbar('eastoutside');
    caxis([zMin, zMax]);
    cbar.Label.String = 'Roughness';
end

%=======================3==================================================
%                                 stack PSDs of roughness profiles
%==========================================================================
[Ny,Nx] = size( ZZ);
% rms roughness
sumZ = sum(sum( ZZ_detr.^2));
rmsZ =  ( 1./(Nx*Ny)*sumZ)^0.5;
sprintf( "rms-roughness: %.3f", rmsZ)
% PSD
m_px=zeros(Ny,Nx);
m_py=zeros(Nx,Ny);

freq =(0:Ny-1)./(Ny*dy);
% a_f=(0:1:nfft-1)./(nfft*delta);
% figure(3)
for i=1:Ny-1
    z1=ZZ_detr(:,i);

    % detrend and demean each profile
    a_z2=detrend(z1-mean(z1));
    %FFTz=fft(z2).*dy;
    a_fftz = 2*abs( fft(a_z2))/Ny;
    m_py(i,:)= a_fftz.^2*dy;
    
    % plot detrended roughness profiles
%     plot( XX(i,:), a_z2, 'k-', 'LineWidth', 2)
%     ylim( zMin, zMax)
 

end

freqX =(0:Nx-1)./(Nx*dx);
for i=1:Nx-1
    z1=ZZ_detr(i,:);
    % detrend and demean each profile
    a_z2=detrend(z1-mean(z1));
    %FFTz=fft(z2).*dy;
    a_fftz = 2*abs( fft(a_z2))/Nx;
    m_px(i,:)=abs(a_fftz.^2).*dx;
end

a_Px=mean(m_px);
a_Py=mean(m_py);

figure(3)
clf()
%plot ave. power spectrum over wavelength
loglog(1./freq,a_Py, 'r', 'LineWidth', 3);
hold on
loglog(1./freqX,a_Px, 'b', 'LineWidth', 3);
xlim(gca, [10, 5*1e3]);
legend('Slip Perp.', 'Slip Parall.')
ylabel('PSD (mu^3)')
xlabel('Wavelength (mu)')

disp( 'Select Roughness Exp. Fitting Range')
[x,y] = ginput(2);

%=======================4==================================================
%                  fit roughness exponent
%==========================================================================
%-log-transform and lsq fit
[par, R] = polyfit(  log10(x), log10(y), 1);
gamma = par(1);
%gamma = 1 + 2H
Hurst = (gamma(1)-1)*.5;
sprintf(  "roughness exponent: %.1f, Hurst=%.2f", gamma, Hurst)







