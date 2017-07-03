function BSCCurve = writeBscFile(fname, maxFreq, diam,sosm,soss,sosshear,rhom,rhos,maxang,maxn)

%a function which calculates scattering cross sections for spheres
%according to Faran's theory; the output is a (1+maxang) by (size of x3) matrix
%where rows correspond to angles (0 to maxang) and columns correspond to different values
%of ka (step size and max value of ka can be changed by manipulating the
%initialization of x3); note that the output is the absolute value of the
%scattering length (i.e. square root of the differential scattering cross
%section) times the wavenumber in the host medium (output is therefore
%unitless)

%--example
% writeBscFile('bsc.dat', 40, 45,1540,5570,3374.7,1020,2540,180,50)
% fname = 'bsc.dat';
% maxFreq = 40 * 10^6; % Hz
% diam = 25;%45;
% sosm = 1540;
% soss = 5570;
% sosshear = 3374.7;%3375;
% rhom = 1020;
% rhos = 2540;
% maxang = 180;
% maxn = 50;

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

%frequencies are in units of Hz and diam is units of microns

freqSpacing = 0.01 * 10^6;
freq = [0:freqSpacing:maxFreq];
sos = 1540;
a = diam/2 * 10^-6; %diameter to radius and microns to meters
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
    jx1 = sphbess(n,x1);
    jx2 = sphbess(n,x2);
    jx3 = sphbess(n,x3);
    
    %spherical neumann functions
    nx3 = sphneumm(n,x3);
    
    %derivatives of spherical bessel and neumann functions
    jpx3 = n/(2*n+1)*sphbess(n-1,x3) - (n+1)/(2*n+1)*sphbess(n+1,x3);
    npx3 = n/(2*n+1)*sphneumm(n-1,x3) - (n+1)/(2*n+1)*sphneumm(n+1,x3);
    jpx1 = n/(2*n+1)*sphbess(n-1,x1) - (n+1)/(2*n+1)*sphbess(n+1,x1);
    jpx2 = n/(2*n+1)*sphbess(n-1,x2) - (n+1)/(2*n+1)*sphbess(n+1,x2);
    
    %Implementation of the equations above Equation (29) in Faran Jr JJ. Sound scattering by solid cylinders and spheres. JASA, 1951, 23(4): 405-418.
    tandeltax3 = -jx3./nx3;
    tanalphax3 = -x3.*jpx3./jx3;
    tanbetax3 = -x3.*npx3./nx3;
    tanalphax1 = -x1.*jpx1./jx1;
    tanalphax2 = -x2.*jpx2./jx2;
    
    %Implementation of Equation (30) in Faran Jr JJ. Sound scattering by solid cylinders and spheres. JASA, 1951, 23(4): 405-418.
    term1 = tanalphax1./(tanalphax1+1);
    term2 = (n^2+n)./(n^2+n-1-0.5*x2.^2+tanalphax2);
    term3 = (n^2+n-0.5*x2.^2+2*tanalphax1)./(tanalphax1+1);
    term4 = ((n^2+n)*(tanalphax2+1))./(n^2+n-1-0.5*x2.^2+tanalphax2);
    tanxsi = -x2.^2/2.*(term1-term2)./(term3-term4);
    
    %Implementation of the equation below Equation (28) in Faran Jr JJ. Sound scattering by solid cylinders and spheres. JASA, 1951, 23(4): 405-418.
    taneta = tandeltax3.*(-rhom/rhos*tanxsi + tanalphax3)./(-rhom/rhos*tanxsi + tanbetax3);
    
    %Implementation of "(2n+1) * sineta * exp(j*eta)" in Equation (31) in Faran Jr JJ. Sound scattering by solid cylinders and spheres. JASA, 1951, 23(4): 405-418. 
    %Note: taneta./(1+taneta.^2) = sineta * coseta
    %      taneta.^2./(1+taneta.^2) = sineta^2
    coeffres = (2*n+1)*taneta./(1+taneta.^2) + i*(2*n+1)*taneta.^2./(1+taneta.^2);
    coeff(n+1,:) = coeffres;
    
    %Implementation of "Pn(costheta)" in Equation (31) in Faran Jr JJ. Sound scattering by solid cylinders and spheres. JASA, 1951, 23(4): 405-418.
    %legendre polynomials
    temp = legendre(n, cos(pi/180*theta));
    
    legend(:,n+1) = temp(1,:);
end

%Refer to Equation (31) in Faran Jr JJ. Sound scattering by solid cylinders and spheres. JASA, 1951, 23(4): 405-418.
k3absscatlength = abs(legend*coeff); %matrix mult completes summation over n

disp('THe number of NANs is: ');
disp(nnz(isnan(k3absscatlength)))
k3absscatlength( isnan(k3absscatlength) ) = 0;

%Absolute scattering length to backscatter coefficients
%Refer to Equation (31) in Faran Jr JJ. Sound scattering by solid cylinders and spheres. JASA, 1951, 23(4): 405-418.
BSCCurve = k3absscatlength./k; %(k3absscatlength./k).^2;

%Write to file
h = fopen(fname, 'wb');
fwrite(h, freqSpacing , 'double');
fwrite(h, length(BSCCurve), 'double' );
fwrite(h, BSCCurve, 'double');
fclose(h);

function sphbess = sphbess(order,vect)

%calculates the spherical bessel function of order 'order' for the passed
%vector

sphbess=sqrt(pi/2./vect).*besselj(order+0.5,vect);

function sphneumm = sphneumm(order,vect)

%calculates the spherical neumann function of order 'order' for the passed
%vector

sphneumm=sqrt(pi/2./vect).*bessely(order+0.5,vect);
