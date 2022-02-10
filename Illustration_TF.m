close all;
clear all;
N_sensors=128;
n=0:127;

%addpath('D:\D\win64_bin\win64_bin');
addpath('D:\tfsa_5-5\windows\win64_bin');
%addpath('E:\Published Papers\DOA estimation of intersecting components 2018\Matlab code');
%addpath('E:\Published Papers\DOA ESTIMATION VITERBI\Multi-sensor IF estimation code');

%crossing componentsi8

s1=exp(2*pi*1i*(0.05*n+0.2*n.^3/(128*128*3)));
s2=1*exp(2*pi*1i*(0.4*n-0*0.3*n.^3/(128*128*3)));
s3=1*exp(2*pi*1i*(0.1*n-0*0.45*n.^3/(128*128*3)));
s4=1*exp(2*pi*1i*(0.475*n-0.2*n.^3/(128*128*3)));
%s5=1*exp(2*pi*1i*(0.35*n-1*n.^2/(128*8)));
%s6=1*exp(2*pi*1i*(0.25*n+0.5*n.^2/(128*8)));
s5=1*exp(2*pi*1i*(0.25*n+0*0.5*n.^2/(128*8)));
%s8=1*exp(2*pi*1i*(0.2*n-1*0.5*n.^2/(128*8)));

perc=0.4;

IF_O(1,:)=0.05+0.2*3*n.^2/(128*128*3);
IF_O(2,:)=0.4-0*0.2*3*n.^2/(128*128*3);
IF_O(3,:)=0.1-0*0.45*3*n.^2/(128*128*3);
IF_O(4,:)=0.475-1*0.2*3*n.^2/(128*128*3);
%IF_O(5,:)=0.35-2*n/(128*8);
%IF_O(6,:)=0.25+2*0.5*n/(128*8);
IF_O(5,:)=0.25+0*0.5*n/(128*8);
%IF_O(8,:)=0.2-2*0.5*n/(128*8);

IF_O=IF_O.';
s = [(s1.') (s2.') (s3.') (s4.') (s5.')];%  (s5.') (s6.') (s7.') ];


clear IF_O;
s1=exp(2*pi*1i*(0.1*n+1*0.2*n.^2/(128*2)));
s2=1*exp(2*pi*1i*(0.3*n-1*0.2*n.^2/(128*2)));
s3=1*exp(2*pi*1i*(0.1*n+1*0.45*n.^3/(128*128*3)));
s4=1*exp(2*pi*1i*(0.475*n-0.2*n.^3/(128*128*3)));
s5=1*exp(2*pi*1i*(0.25*n+1*0.75*n.^2/(128*8)));

perc=0.4;

IF_O(1,:)=0.05+0.2*3*n.^2/(128*128*3);
IF_O(2,:)=0.4-1*0.3*3*n.^2/(128*128*3);
IF_O(3,:)=0.1+1*0.45*3*n.^2/(128*128*3);
IF_O(4,:)=0.475-1*0.2*3*n.^2/(128*128*3);
IF_O(5,:)=0.25+1.5*n/(128*8);
%IF_O(6,:)=0.25+2*0.5*n/(128*8);
%IF_O(8,:)=0.2-2*0.5*n/(128*8);

IF_O=IF_O.';
s = [(s1.') (s2.') (s3.') (s4.') (s5.')];%  (s5.') (s6.') (s7.') ];
s = [(s1.') (s2.') ];

n_sources=2;



s_orig=s;

% set mixing matrix A
%theta = [15,30,50]*pi/180;   % sensor separation angles in radians
%        theta = [-12,0,12]*pi/180;   % sensor separation angles in radians
theta = [0,10,20,30,40,50,60]*pi/180;   % sensor separation angles in radians
theta = [0,10,20,30,40]*pi/180;   % sensor separation angles in radians
theta = [0,30]*pi/180;   % sensor separation angles in radians

LL=100;
index=0;
delta=2;
num=n_sources;
        
        A = exp(1j*pi*[0:N_sensors-1].'*sin(theta));  % mixing matrix A
        
        
        X = A*s.';                             % mixed source
        theta9=round(theta *180/pi);
        % generate noise
        SNR=100;
        sigma = 10^(-SNR/20);
        w = sigma*(randn(N_sensors,length(n)) + 1j*(randn(N_sensors,length(n))))/sqrt(2); % noise
        
        X=X+w;
        
for i=1:N_sensors
   I(i,:,:)=imresize(spec(X(i,:),1,35,'hamm',128,1),[128,128]); 
end

figure;
imagesc(reshape(abs(I(1,:,:)),128,128))
set(gcf,'Position',[20 100 640 500]);
xlabel('Time / Sec','FontSize',20,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'YDir','normal');
set(gca,'FontSize',20);
[x,y,z]=meshgrid(1:64,1:64,1:64);
II=abs(fft(I,[],1));
%II=I;
idx = find(II/max(II(:))>=0.6);
[X, Y, Z] = ind2sub(size(I), idx);
pointsize = 30;
%scatter3(X(:), Y(:), Z(:), pointsize, A(idx));
figure;
scatter3( X(:), Y(:),Z(:), pointsize, II(idx));
xlabel('Spatial Frequency');
ylabel('Time');
zlabel('Frequency');



%imagesc(reshape(abs(II(64,:,:)),64,64))
