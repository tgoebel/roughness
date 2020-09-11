% 
% fault roughness analysis, based on 2D surface scans
%
% 
%
% @author thgoebel, CERI U of Memphis

%=======================1==================================================
%                  parameters
%==========================================================================
data_dir = "~/Documents/teaching/CERI/DataAnalysis_7104/data/";
file_in  = "fault_roughness.txt";

zMin = -40;
zMax =  40;

XYZ = load(strcat(data_dir,file_in));
dx  = XYZ(2,1) - XYZ(1,1);
fprintf( '%s%.2f', 'grid spacing :', dx)
%=======================2==================================================
%                  data pre-processing
%==========================================================================
%A% remove nans
XYZ =  XYZ(~isnan(XYZ(:,3)),:);
Z   =  XYZ(:,3);
X   =  XYZ(:,1);
Y   =  XYZ(:,2);
sprintf( '%s %.1f %.1f','range of topography:',min(Z), max(Z))
%B% interpolation, fill gaps
a_x = unique( X);
a_y = unique( Y);
[XX, YY] = meshgrid( a_x, a_y);
size(XX)
ZZ       = griddata(X, Y, Z, XX, YY, 'linear');

%C%visualize raw data
figure(1)
clf
pcolor( XX, YY, ZZ); 
shading interp;
cbar = colorbar('eastoutside');
cbar.Label.String = 'Roughness';

%D% detrend
ZZ_detr = detrend_2d(ZZ);
figure(2)
clf
pcolor( XX, YY, ZZ_detr); 
colormap hot
shading flat;
cbar = colorbar('eastoutside');
caxis([zMin, zMax]);
cbar.Label.String = 'Roughness';

%=======================3==================================================
%                  stack roughness profiles
%==========================================================================
[Nx,Ny] = size( ZZ);
m_px=zeros(Ny,Nx);
m_py=zeros(Nx,Ny);


freq =(0:Nx-1)./(Nx*dx);
for i=1:Nx-1
    z1=ZZ_detr(:,i)';
    if ~sum(isnan(z1))
        mean( z1)
        z1=z1-mean(z1);
        z2=detrend(z1);
        FFTz=fft(z2).*dx;
        m_px(i,:)=abs(FFTz.^2).*(Nx*dx);
%         loglog(1./freq,m_px(i,:));
%         xlim(gca, [10, 5*1e3]); 
%         hold on
    else
        m_px(i,:)=nan;
    end
end
a_Py=nanmean(m_px);

figure(3)
%plot ave. power spectrum over wavelength
loglog(1./freq,a_Py, 'r', 'LineWidth', 3);
xlim(gca, [10, 5*1e3]);
legend('Perp')
ylabel('PSD (mu^3)')
xlabel('Wavelength (mu)')

disp( 'Select Roughness Exp. Fitting Range')
[x,y] = ginput(2);

%=======================4==================================================
%                  fit roughness exponent
%==========================================================================
%TODO: -log-transform and lsq fit

gamma = log10(x)\log10(y);
sprintf( '%s%.1f%s', "roughness exponent: ", gamma)

%TODO: - determine difference in fitting exponent for different bandwidth (wavelength)





