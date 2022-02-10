

function [fidexmult] = FASTEST_IF(Sig,win_length, num, delta,L,thr,Thr,step,FFT_len)
% SET OF FRACTIONAL WINDOWS
w=gausswin(win_length,1);
%w=ones(win_length,1);
l=0:win_length-1;
%FFT_len=win_length-1;
%if nargin==7
%step=round(win_length/4);
%end
for iii=1:FFT_len
    WW(iii,:)=exp(-1i*(iii)*2*pi*l/FFT_len);
end
e_max=0;
%tic;
i=0;
window_rot=zeros(2*L+1,win_length);
for k=-L+1:1:L-1
   i=i+1;
   window_rot(i,:)=frft(w,1*k/L);%fracft(w,0.95* k/L);%0.05
%   window_rot(i,:)=fracft(w,k/L);%0.05
end
%save('window_rot','window_rot');
%load('window_rot');


w_signal=zeros(1,win_length);
v=zeros(1,2*L+1);
index=v;

for iii=1:num
%I=quadtfd(Sig,length(Sig)-1,1,'specx',35,'hamm',256);
    Sig_extended=[zeros(1,floor(win_length/2)) Sig zeros(1,floor(win_length/2))];

    Siga=filter(ones(1,win_length),1,abs(Sig));
    [~,t_start]=max(Siga(floor(win_length/2)+1:end-floor(win_length/2)));
    t_start=t_start(1)+floor(win_length/2);

    for i=1:2*L+1
        FF=abs(fft(Sig(t_start-floor(win_length/2):t_start+floor(win_length/2)).*window_rot(i,:),FFT_len));
        
        [v(i),index(i)]=max(FF(1:round(end/2)));
      %          [v(i),index(i)]=max(FF(1:end));

    end
    [v_m,ind]=max(v);
    
    freq_start=index(ind)-1;
    frac_wind_index_start=ind;
    v_oldd=v_m;
    
    
    IF=zeros(1,length(Sig));
    IF(t_start)=freq_start;
    
    
    clear v;
    

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

            
           k=f0-delta:1:f0+delta;
           k(k>FFT_len/2)=FFT_len/2;
           k(k<=0)=1;
            if frac_wind_index<2
                frac_wind_index=2;
            elseif frac_wind_index>2*L-1
                
                frac_wind_index=2*L-1;
            end
            %   for i=frac_wind_index-1:frac_wind_index+1    % FOR ALL WINDOWS
            w_signal(1:win_length)=(Sig_extended(t_index:t_index+win_length-1).*window_rot(frac_wind_index-1,:)); %WINDOWED SIGNAL
            V=abs(WW(k,:)*w_signal.');
            [v(1),indexx(1)]=max(V);
            
            w_signal(1:win_length)=(Sig_extended(t_index:t_index+win_length-1).*window_rot(frac_wind_index,:)); %WINDOWED SIGNAL
            V=abs(WW(k,:)*w_signal.');
            [v(2),indexx(2)]=max(V);
            
            w_signal(1:win_length)=(Sig_extended(t_index:t_index+win_length-1).*window_rot(frac_wind_index+1,:)); %WINDOWED SIGNAL
            V=abs(WW(k,:)*w_signal.');
            [v(3),indexx(3)]=max(V);
            [v_m,ind]=max(v);
            f0=k(indexx(ind));%+f0-delta;
            frac_wind_index=frac_wind_index+ind-2;
            
            
            if v_m<Thr*v_oldd
                
                break;
            end
            %v_old=v_m;
            
            IF(t_index)=f0;
        end
    end
        ind=find(IF>0);
        IF=interp1(ind,IF(ind),1:length(Sig));
      IF(isnan(IF))=0;
    IF=IF/(1*FFT_len);
    Phase=2*pi*filter(1,[1 -1],IF);
    s_dechirp=exp(-1i*Phase);
    
    
    LL=delta;
    %TF filtering for each sensor
    s1 = Sig.*(s_dechirp);
    s2=fftshift(fft(s1));
    %figure; plot(abs(s2));
    if sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2))<e_max*thr
        break;
    else
        %  e_max
        if sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2))>e_max
            e_max=sum(abs(s2(length(Sig)/2-LL-1:length(Sig)/2+LL-1).^2));
        end
        % Energy of the last component
        s2(length(Sig)/2-LL:length(Sig)/2+LL)=0;
        s2=ifft(ifftshift(s2)).*conj(s_dechirp);
        Sig=s2;%-extr_Sig(iii);
        fidexmult(iii,:) = IF;
    end
    
end


end


