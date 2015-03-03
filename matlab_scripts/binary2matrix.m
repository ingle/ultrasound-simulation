function A0 = binary2matrix(rf_file)
%%%%%the RF data is read from file and converted to B-mode values
%%Adjust the center frequency and bandwidth%%%
bandwidth=.5;
centerFreq=5e6;   %in Hz

[A0, fs]=generate_rf( rf_file, centerFreq, bandwidth);

%info for plotting
BeamSpace=0.2; %again in mm, equal to element spacing?
zDim= 42;  %of beam, in mm
nLine = 100;
pixelThick=zDim/(size(A0,1)-1);


x=((1:nLine)-(nLine+1)/2)*BeamSpace;
y=( 1:size(A0,1) )*pixelThick;
logenv=20*log10(abs(A0)+1e-5);
logenv=logenv-max(logenv(:));


% %the B-mode data is plotted
imagesc(x, y, logenv, [-60, 0]); axis equal;  axis image,axis tight;colormap(gray);
A0 = real(A0);





function [rfdata, fs]=generate_rf(fname, centerfreq, bw)
%%Input:
%fname = name of the file containing rf data output by simulation
%centerfreq = Desired center frequency of the simulated transducer in units of Hz
%bw = The bandwidth of the transducer

%%open and read simulation file
fp=fopen(fname,'r');
%first read in the deltaF, lines, and points
freqstep = fread(fp, 1, 'double');
freqpoints = fread(fp, 1, 'int');
nLines = fread(fp, 1, 'int' );

fs = freqstep*freqpoints;

%%%read in real and imaginary parts of RF data in frequency domain
tempReal=fread(fp, [(freqpoints), nLines],'double');
tempImag =fread(fp, [(freqpoints), nLines],'double');
freqdata=tempReal - 1i*tempImag;
fclose(fp);
clear tempReal tempImag;


%%%create transducer spectrum
alpha=4*log(2)/(centerfreq*centerfreq*bw*bw);

f=(0:freqpoints)*freqstep;
pulsespectrum=1e6*exp(-alpha*(f-centerfreq).*(f-centerfreq) ) ;

[row, col]=size(freqdata);
rfSpectrum=freqdata.*(ones(col,1)*pulsespectrum(1:row))';

%zero-pad the data prior to transform to achieve higher sampling frequency
rfdata=ifft(rfSpectrum, freqpoints*3);
