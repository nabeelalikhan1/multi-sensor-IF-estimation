

function [fidexmult,Xout] = Multi_Sensor_FASTEST_IF_new1(Sig,N_S,win_length, num, delta,L,thr,Thr,step,FFT_len)
% SET OF FRACTIONAL WINDOWS
% Sig=X;
% N_S=N_sensors;
% win_length=65;
% num=n_sources;
% thr=0;
% Thr=0;
% L=8;
% FFT_len=length(X);
% step=4;
%delta
w=gausswin(win_length,1);
l=0:win_length-1;
FFT_len1=max(2*N_S,64);
kl=0:N_S-1;
for iii=1:FFT_len1
    W2(iii,:)=exp(-1i*(iii)*2*pi*kl/FFT_len1);
end

for iii=1:FFT_len
    WW(iii,:)=exp(-1i*(iii)*2*pi*l/FFT_len);
end
e_max=0;
%tic;
i=0;
window_rot=zeros(2*L+1,win_length);
for k=-L+1:1:L-1
    i=i+1;
    window_rot(i,:)=frft(w,1* k/L);%0.05
end
%save('window_rot','window_rot');
%load('window_rot');
%toc;

for iii=1:num
    for jj=1:N_S
        Sig_extended(jj,:)=[zeros(1,floor(win_length/2)) Sig(jj,:) zeros(1,floor(win_length/2))];
    end
    Siga=filter(ones(1,win_length),1,sum(abs(Sig)));
    % Siga=
    [~,t_start]=max(Siga(floor(win_length/2)+1:end-floor(win_length/2)));
    t_start=t_start(1)+floor(win_length/2);
    % t_start=46;
    v_m=0;
    for i=1:2*L+1
        FF=zeros(N_S,length(Sig));
        for jj=1:N_S
            A=(fft(Sig(jj,t_start-floor(win_length/2):t_start+floor(win_length/2)).*window_rot(i,:),1*length(Sig)));
            FF(jj,1:end/2)=A(1:end/2);
        end
        FFd=abs(fft(FF,FFT_len1));
        v=max(FFd(:));
        [indexS,index]=find(FFd==max(FFd(:)));
        %index=f(1);
        %f1=f1(1);
            %FF=abs(fft(fft(Sig(:,t_start-floor(win_length/2):t_start+floor(win_length/2)).*window_rot(i,:),length(Sig)));
        
        %[v,index]=max(FF(1:end/2));
        if v_m<v
            freq_start=index-1;  %peak_frequency location
            frac_wind_index_start=i;   %chirp rate location
            freqS=indexS;
            v_m=v;
        end
    end
    v_oldd=v_m;
    
    %t_start=t_start-floor(win_length/2)-1;
    IF=zeros(1,length(Sig))-1;
    IF(t_start)=freq_start;
    %freqS
   % frac_wind_index_start=L;
    for it=1:2
        f0=freq_start;
        frac_wind_index=frac_wind_index_start;
        t_index=t_start;
        while and(t_index>1,t_index<length(Sig))
            if it==1
                
                if t_index+step>length(Sig)
                    %f0=f0+(length(Sig)/win_length)*(length(Sig)-t_index)*tan(1*pi*(frac_wind_index-(L))/(2*L));
                    f0=f0+(FFT_len/win_length)*(length(Sig)-t_index)*tan(1*pi*(frac_wind_index-(L))/(2*L));
                    t_index=length(Sig);
                    %             IF(t_index)=round(f0);
                    
                    %  break;
                else
                    t_index=t_index+step;
                    f0=f0+(FFT_len/win_length)*step*tan(1*pi*(frac_wind_index-(L))/(2*L));
                end
                
            else
                if t_index-step<1
                    f0=f0-(FFT_len/win_length)*abs(1-t_index)*tan(1*pi*(frac_wind_index-(L))/(2*L));
                    t_index=1;
                    %            IF(t_index)=round(f0);
                    %break;
                    
                else
                    t_index=t_index-step;
                    f0=f0-(FFT_len/win_length)*step*tan(1*pi*(frac_wind_index-(L))/(2*L));
                end
                
                
            end
            f0=round(f0);
            
            v_m=0;
            k=f0-delta:1:f0+delta;
            k(k>FFT_len/2)=FFT_len/2;
            k(k<=0)=1;
            kj=freqS;
            kj(kj<=0)=1;
            kj(kj>FFT_len1)=FFT_len1;
            %kj
            %kj=indexS;
            %k33=freqS
            if frac_wind_index<2
                frac_wind_index=2;
            elseif frac_wind_index>2*L-1
                
                frac_wind_index=2*L-1;
            end
            for i=frac_wind_index-1:frac_wind_index+1    % FOR ALL WINDOWS

            for jj=1:N_S
                w_signal(jj,1:win_length)=(Sig_extended(jj,t_index:t_index+win_length-1).*window_rot(i,:)); %WINDOWED SIGNAL
            end
            for jjj1=1:length(kj)
            for jjj=1:length(k)%jj
                VVV=0;
                for jj=1:N_S
                    VVV(jj)=((sum(w_signal(jj,:).*WW(k(jjj),:))));
                end
                
                %V(jjj)=abs(sum(VVV.*W2(freqS,:)));
                V(jjj1,jjj)=abs(sum(VVV.*W2(kj(jjj1),:)));
            end
            end
            v=max(V(:));
           % [v,index]=max(V);
            [f1,index]=find(V==v);
            %[index,f1]=find(V==v);
            
            index=index(1);
            f1=f1(1);
            if v_m<v
                f0=k(index);  %peak_frequency location
                frac_wind_index=i;   %chirp rate location
                v_m=v;
                freqSS=kj(f1);
            end
        end
        if v_m<Thr*v_oldd
            
            break;
        end
        %freqSS
        IF(t_index)=f0;
       % freqS=round(0.8*freqS+0.2*freqSS);
    end
    end
    %IF(t_index)=f0;
    
ind=find(IF>0);
        IF=interp1(ind,IF(ind),1:length(Sig));
      IF(isnan(IF))=0;
    IF=IF/(1*FFT_len);
    
    

%figure; plot(IF)
Phase=2*pi*filter(1,[1 -1],IF);
s_dechirp=exp(-1i*Phase);


LL=delta;
%TF filtering for each sensor
fidexmult(iii,:) = IF;
ssss=0;
for jj=1:N_S
    
    s1 = Sig(jj,:).*(s_dechirp);
    s2=fftshift(fft(s1));
    %figure; plot(abs(s2));
    if sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2))<e_max*thr
        break;
    end
    %  e_max
    if sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2))>e_max
        e_max=sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2));
    end
    ssss=ssss+sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2));
    % Energy of the last component
    s1=zeros(size(s2));
    s1(length(Sig)/2-LL:length(Sig)/2+LL)=s2(length(Sig)/2-LL:length(Sig)/2+LL);
    s2(length(Sig)/2-LL:length(Sig)/2+LL)=0;
    
    s1=ifft(ifftshift(s1)).*conj(s_dechirp);
    s2=ifft(ifftshift(s2)).*conj(s_dechirp);
    Sig(jj,:)=s2;%-extr_Sig(iii);
    Xout(jj,iii,:)=s1;
    % end
end
if ssss<e_max*thr
    break;
end
if ssss>e_max
    e_max=ssss;
end

end



