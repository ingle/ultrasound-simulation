function BSCCurve = sphere(fname, maxFreqMHz, diam,sosm,soss,sosshear,rhom,rhos,maxang,maxn)

%This function writes a backscatter coefficient file given a filename
%and a maximum frequency.


%initialization of x3); note that the output is the absolute value of the
%scattering length (i.e. square root of the differential scattering cross
%section) times the wavenumber in the host medium (output is therefore
%unitless)

%sosm = speed of sound in host medium (fluid)
%soss = speed of sound in sphere
%sosshear = speed of sound for shear waves in sphere (if poisson's ratio
%(sigma)is available, sosshear=soss*sqrt((1-2*sigma)/2/(1-sigma)))
%rhom = host medium density
%rhos = sphere density
%maxang = maximum angle (in degrees) for which to calculate a cross
%section
%maxn = number of terms to include in the summation (higher number =
%greater accuracy but more computation time)

%UNITS ARE OF NO CONSEQUENCE AS LONG AS THEY ARE CONSISTENT, I.E.
%SPEEDS OF SOUND MUST HAVE THE SAME UNITS AND DENSITIES MUST HAVE THE SAME
%UNITS

%9/20/04 Tony Gerig

%Hairong Shi added comments:
% for glass beads we used in lab, Poisson's ratio sigma=0.21,
%sosm=1540m/s
%soss=5570m/s
%so sosshear=3374.7m/s
%rhom=1020 kg/m^3
%rhos=2540 kg/m^3
%
%so an example to run this code is:
% BSCCurve = sphere(freq, diam,1540,5570,3374.7,1020,2540,180,50);

%frequencies are in units of MHz and diam is units of microns

freqSpacingMHz = .01;
freq = [0:freqSpacingMHz:maxFreqMHz];
sos = 1540;
a = diam/2 *10^-6; %diameter to radius and microns to meters
k = 2*pi./(sos./freq);



theta=maxang;
x3=k*a; %the range of ka over which cross sections are calculated
x2=x3*sosm/sosshear; %ka value for shear in sphere
x1=x3*sosm/soss; %ka value for sphere

for n=0:maxn %main loop over order n
    if n == 90
        q = 0;
    end
    %spherical bessel functions
    jx1=sphbess(n,x1);
    jx2=sphbess(n,x2);
    jx3=sphbess(n,x3);
    
    %spherical neumann functions
    nx3=sphneumm(n,x3);
    
    %derivatives of spherical bessel and neumann functions
    jpx3=n/(2*n+1)*sphbess(n-1,x3)-(n+1)/(2*n+1)*sphbess(n+1,x3);
    npx3=n/(2*n+1)*sphneumm(n-1,x3)-(n+1)/(2*n+1)*sphneumm(n+1,x3);
    jpx1=n/(2*n+1)*sphbess(n-1,x1)-(n+1)/(2*n+1)*sphbess(n+1,x1);
    jpx2=n/(2*n+1)*sphbess(n-1,x2)-(n+1)/(2*n+1)*sphbess(n+1,x2);
    
    tandeltax3=-jx3./nx3;
    tanalphax3=-x3.*jpx3./jx3;
    tanbetax3=-x3.*npx3./nx3;
    tanalphax1=-x1.*jpx1./jx1;
    tanalphax2=-x2.*jpx2./jx2;
    
    %calculation of tanxsi and coefficient
    term1=tanalphax1./(tanalphax1+1);
    term2=(n^2+n)./(n^2+n-1-0.5*x2.^2+tanalphax2);
    term3=(n^2+n-0.5*x2.^2+2*tanalphax1)./(tanalphax1+1);
    term4=((n^2+n)*(tanalphax2+1))./(n^2+n-1-0.5*x2.^2+tanalphax2);
    tanxsi=-x2.^2/2.*(term1-term2)./(term3-term4);
    
    taneta=tandeltax3.*(-rhom/rhos*tanxsi+tanalphax3)./(-rhom/rhos*tanxsi+tanbetax3);
    coeff(n+1,:)=(2*n+1)*taneta./(1+taneta.^2)+i*(2*n+1)*taneta.^2./(1+taneta.^2);
    
    %legendre polynomials
    temp=legendre(n,cos(pi/180*theta));
    legend(:,n+1)=temp(1,:);
end

k3absscatlength=abs(legend*coeff); %matrix mult completes summation over n

%absolute scattering length to backscatter coefficients
disp('THe number of NANs is: ');
disp(nnz(isnan(k3absscatlength)))
k3absscatlength( isnan(k3absscatlength) ) = 0;
BSCCurve=(k3absscatlength./k).^2;


h = fopen(fname, 'wb');
fwrite(h, freqSpacingMHz , 'double');
fwrite(h, length(BSCCurve), 'double' );
fwrite(h, BSCCurve, 'double');

function sphbess = sphbess(order,vect)

%calculates the spherical bessel function of order 'order' for the passed
%vector

sphbess=sqrt(pi/2./vect).*besselj(order+0.5,vect);

function sphneumm = sphneumm(order,vect)

%calculates the spherical neumann function of order 'order' for the passed
%vector

sphneumm=sqrt(pi/2./vect).*bessely(order+0.5,vect);
