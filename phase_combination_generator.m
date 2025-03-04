%This code is to check the effects of the third order correction
clc;clear;close all;
etah=figure;
poth=figure;
L=4097;
% L=8193;
Tp=12;
e=exp(1);
wp=2*pi/Tp;
RA_on=1;
ORDER=5;
% phasecomb=90:5:270;
% phases=[0,phasecomb];
phasecomb1=5:5:85;
phasecomb2=275:5:355;
phases=[phasecomb1,phasecomb2];
% Alpha1=7:-1:2;
% Alpha2=1.9:-0.1:1;
% Alpha_All=[Alpha1, Alpha2];
Alpha_All=[7];
% for ORDER=1:5
for Akp=0.02:0.02:0.18
    for Alpha=Alpha_All
        for phase_index=1:length(phases)
            phi_shift=phases(phase_index);
            %         write_path=fullfile(pwd,"test1",sprintf("gamma_%d_Akp_%.3d_phi_%d",gamma,round(Akp*100),phi_shift));
            %             write_path=fullfile(pwd,"test7",sprintf("Alpha_%.1f_Akp_%.3d_phi_%d",Alpha,round(Akp*100),phi_shift));

            t=-40*Tp;
            %         h=54;%Shallow Water
            h=150;
            %Deep water
            %         h=32;%kd=0.88
            %         h=54;%kd=1.5
            g=9.81;
            % kp = abs(solveWForK(wp, g, h));
            kp=wp.^2/g;
            kd=kp*h;
            lambda=2*pi/kp;

            A=Akp/kp;
            Nx=68*lambda;
            %             dx=lambda/30;
            dx=lambda/30;
            N=Nx/dx;
            kmax=2*pi/dx/2;

            dk=kmax/N;
            k=dk:dk:kmax;
            T_total=60*Tp;
            dt=Tp/30;
            N_steps=T_total/dt;

            cp=sqrt(g/kp*tanh(kp*h));
            cw=1/2*cp*(1+2*kp*h/sinh(2*kp*h));
            CFL=cw*dt/dx;
            %     %%%%%%%%%%%%%%%%%%%Using Gaussian Spectrum
            %             kw=BandM*0.004606;
            %
            %             for i=1:numel(k)
            %                 S(i)=exp(-(k(i)-kp).^2/(2*kw.^2));
            %             end
            %             %         semilogy(k/kp,S);
            %             %         ylim([1e-6,1e1]);
            %             for i=1:numel(k)
            %                 a(i)=A*S(i)*dk/(sum(S*dk));
            %             end
            %             % a=sqrt(2*S*dk*2*pi);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Using Semi-Gaussian
            kw=0.004606;
            for i=1:numel(k)
                if k(i)<kp
                    S(i)=exp(-(k(i)-kp).^2/(2*kw.^2));
                else
                    %        S(i)=exp(-(k(i)-kp).^2/(2*kw.^2))^(1/alpha);
                    kw2=sqrt(kp^2/(2*log(10^Alpha)));
                    S(i)=exp(-(k(i)-kp).^2/(2*kw2.^2));
                end
            end
            %         semilogy(k/kp,S);
            %         ylim([1e-6,1e1]);
            for i=1:numel(k)
                a(i)=A*S(i)*dk/(sum(S*dk));
            end
            w=sqrt(g.*k.*tanh(k.*h));
            %     %%%%%%%%%%%%%%%%%%%%%%%%%%JONSWAP
            %             w=sqrt(g.*k.*tanh(k.*h));
            %             % wp=sqrt(g.*kp.*tanh(kp.*h));
            % %             gamma=3.3;
            %             for i=1:numel(k)
            %                 if k(i)<kp
            %                     sigma=0.07;
            %                 else
            %                     sigma=0.09;
            %                 end
            %                 % sigma=0.2;
            %                 alpha=0.0624/(0.23+0.0336*gamma-0.185/(1.9+gamma));
            %                 S(i)=4*pi^2*alpha/(2*k(i)^3)*exp(-1.25*(kp/k(i))^4)*gamma^(exp(-(k(i)-kp)^2/(2*sigma^2*kp^2)));
            %             end
            %
            %
            %             for i=1:numel(k)
            %                 a(i)=A*S(i)*dk/(sum(S*dk));
            %             end
            %             semilogy(k/kp,a);
            %             ylim([1e-6,1e1]);



            phi=-w*t+k.*t*cw;
            phi=phi+deg2rad(phi_shift);
            % phi=phi';
            %Linear Eta
            % a=a*numel(a);
            % F0=complex(a.*cos(phi)*numel(a),a.*sin(phi)*numel(a));
            F0=a.*numel(a).*exp(1i.*phi);
            F=complex([0,real(F0),fliplr(real(F0))],[0,imag(F0),-fliplr(imag(F0))]);
            XX_linear=ifft(F);
            XX_linear=fftshift(XX_linear);

            F0=a.^2.*numel(a).*(k)./(4.*tanh(k.*h)).*(2+3./(sinh(k.*h).^2));
            F0=a20array(F0);
            phi2=a20array(phi);
            F0=F0.*exp(1i.*phi2);
            F=complex([0,real(F0),fliplr(real(F0))],[0,imag(F0),-fliplr(imag(F0))]);
            XX_20=ifft(F);
            XX_20=fftshift(XX_20);
            XX_20=XX_20-sum(a.^2.*numel(a).*k/(2*sinh(2.*k*h)));
            %Linear Potential
            for i=1:numel(k)
                ap(i)=g/(w(i))*A*S(i)*dk/(sum(S*dk));%Valid for deep water only
            end
            ap=ap*numel(ap);
            % Fp=complex([0,real(Fp0),fliplr(real(Fp0))],[0,imag(Fp0),-fliplr(imag(Fp0))]);
            PP_linear=zeros(1,numel(XX_linear));
            for i=1:numel(XX_linear)
                z=XX_linear(i);
                Fp=ap.*cosh(k.*(z+h))./cosh(k.*h);
                Fp=Fp(1,:);
                Fp0=Fp.*exp(1i*(phi-0.5*pi));
                Fp=complex([0,real(Fp0),fliplr(real(Fp0))],[0,imag(Fp0),-fliplr(imag(Fp0))]);
                PPP=ifft(Fp);
                PPP=fftshift(PPP);
                PP_linear(i)=PPP(i);

                Fp20=a.^2.*numel(a).*3.*w/8.*cosh(2.*k.*(z+h))./(sinh(k*h).^4);
                Fp20=a20array(Fp20);
                Fp20=Fp20.*exp(1i*(phi2-0.5*pi));
                Fpp20=complex([0,real(Fp20),fliplr(real(Fp20))],[0,imag(Fp20),-fliplr(imag(Fp20))]);
                PPP=ifft(Fpp20);
                PPP=fftshift(PPP);
                PP_20(i)=PPP(i);
            end


            % z=0;
            % Fp=ap.*cosh(k.*(z+h))./cosh(k.*h);
            % Fp0=Fp.*exp(1i*(phi-0.5*pi));
            % Fp=complex([0,real(Fp0),fliplr(real(Fp0))],[0,imag(Fp0),-fliplr(imag(Fp0))]);
            % PPP=ifft(Fp);
            % PPP=fftshift(PPP);

            %             if ORDER~=1
            %Second Order Coeffs
            mu1=1;
            mu2=1;
            Bp=zeros(numel(S),numel(S));
            Bm=zeros(numel(S),numel(S));
            Ap=zeros(numel(S),numel(S));
            Am=zeros(numel(S),numel(S));
            Dp=zeros(numel(S),numel(S));
            Dm=zeros(numel(S),numel(S));

            for i=1:numel(S)
                for j=1:numel(S)
                    if i~=j
                        Dp(i,j)=(w(i)+w(j)).^2-g*abs(k(i)+k(j))*tanh(abs(k(i)+k(j))*h);
                        Dm(i,j)=(w(i)-w(j)).^2-g*abs(k(i)-k(j))*tanh(abs(k(i)-k(j))*h);
                        Ap(i,j)=-w(i)*w(j)*(w(i)+w(j))/Dp(i,j)*(1-cos(mu1-mu2)/(tanh(k(i)*h)*tanh(k(j)*h)))+1/(2*Dp(i,j))*(w(i).^3/(sinh(k(i)*h)^2)+w(j).^3/(sinh(k(j)*h)).^2);
                        Am(i,j)=w(i)*w(j)*(w(i)-w(j))/Dm(i,j)*(1+cos(mu1-mu2)/(tanh(k(i)*h)*tanh(k(j)*h)))+1/(2*Dm(i,j))*(w(i).^3/(sinh(k(i)*h)^2)-w(j).^3/(sinh(k(j)*h)).^2);
                        Bp(i,j)=(w(i).^2+w(j).^2)/(2*g)-w(i)*w(j)/(2*g)*(1-cos(mu1-mu2)/(tanh(k(i)*h)*tanh(k(j)*h)))*(((w(i)+w(j))^2+g*abs(k(i)+k(j))*tanh(abs(k(i)+k(j))*h))/(Dp(i,j)))+(w(i)+w(j))/(2*g*Dp(i,j))*(w(i)^3/(sinh(k(i)*h)^2)+w(j)^3/(sinh(k(j)*h)^3));
                        Bm(i,j)=(w(i).^2+w(j).^2)/(2*g)+w(i)*w(j)/(2*g)*(1+cos(mu1-mu2)/(tanh(k(i)*h)*tanh(k(j)*h)))*(((w(i)-w(j))^2+g*abs(k(i)-k(j))*tanh(abs(k(i)-k(j))*h))/(Dm(i,j)))+(w(i)-w(j))/(2*g*Dm(i,j))*(w(i)^3/(sinh(k(i)*h)^2)-w(j)^3/(sinh(k(j)*h)^3));
                    end
                end
            end

            %Second Order Eta
            X2m=zeros(1,numel(XX_linear));
            X2p=zeros(1,numel(XX_linear));
            xx=1:numel(XX_linear);
            xx=xx*dx;
            xx=xx-xx(end)/2;
            for i=1:numel(a)
                %Second Order Sum Correction for Eta
                Bpn=a(i).*a.*Bp(i,:);
                Bpn=Bpn*numel(Bpn);
                Bpn=circshift(Bpn,i);
                Phii=phi+phi(i);
                Phii=circshift(Phii,i);
                Fn=Bpn.*exp(1i*Phii);
                Fn=complex([0,real(Fn),fliplr(real(Fn))],[0,imag(Fn),-fliplr(imag(Fn))]);
                XXn=ifft(Fn);
                XXn=fftshift(XXn);
                X2p=X2p+XXn;
                %Second Order Sub Correction for Eta
                %%% phi=w*t-k.*t*cw;


                %     Bmn=a(i).*a.*Bm(i,:);
                %     Phii=phi-phi(i);
                %     Bmn=circshift(Bmn,-i);
                %     Phii=circshift(Phii,-i);
                Bmnp=zeros(1,numel(Bpn));
                Bmnm=zeros(1,numel(Bpn));
                Phii1=zeros(1,numel(Bpn));
                Phii2=zeros(1,numel(Bpn));
                for j=1:numel(Bmnp)
                    if i-j>0
                        Bmnp(i-j)=a(i)*a(j)*Bm(i,j);
                        Phii1(i-j)=phi(i)-phi(j);
                    end
                    if i-j<0
                        Bmnm(j-i)=a(i)*a(j)*Bm(i,j);
                        Phii2(j-i)=-phi(i)+phi(j);
                    end
                end


                Bmnp=Bmnp*numel(Bmnp);
                % Fn=complex(Amn.*cos(Phii),Amn.*sin(Phii));

                Fn=Bmnp.*exp(1i*(Phii1));
                Fn=complex([0,real(Fn),fliplr(real(Fn))],[0,imag(Fn),-fliplr(imag(Fn))]);
                XXn=ifft(Fn);
                XXn=fftshift(XXn);
                X2m=X2m+XXn;

                Bmnm=Bmnm*numel(Bmnm);
                Fn=Bmnm.*exp(1i*(Phii2));
                Fn=complex([0,real(Fn),fliplr(real(Fn))],[0,imag(Fn),-fliplr(imag(Fn))]);
                XXn=ifft(Fn);
                XXn=fftshift(XXn);
                X2m=X2m+XXn;

            end
            X2p=X2p/2;
            X2m=X2m/2;
            % hold on
            % plot(X2p);
            % plot(X2m);
            %Remove the error from ifft
            N=length(X2p);
            X2p(1:round(0.05*N))=0;
            X2p(round(0.95*N):end)=0;
            N=length(XX_20);
            XX_20(1:round(0.1*N))=0;
            XX_20(round(0.9*N):end)=0;


            % dx=pi/k_max;
            XX_2=XX_linear+X2m+X2p+XX_20;

            %Second Order Potential
            P2m=zeros(1,numel(PP_linear));
            P2p=zeros(1,numel(PP_linear));
            %             for X_index=1:numel(XX_linear)
            %                 z=XX_linear(X_index);
            z=0;
            P2m=zeros(1,numel(PP_linear));
            P2p=zeros(1,numel(PP_linear));
            for i=1:numel(a)
                %Second Order Sum Correction for Potential
                %Assume z=XX_linear(i) and performing ifft several times
                Apn=a(i).*a.*Ap(i,:).*cosh((k(i)+k).*(z+h))./cosh((k(i)+k)*h);
                Apn=Apn*numel(Apn);
                Apn=circshift(Apn,i);
                Phii=phi+phi(i);
                Phii=circshift(Phii,i);


                % Fn=complex(Apn.*cos(Phii),Apn.*sin(Phii));
                Fn=Apn.*exp(1i*(Phii-0.25*pi));
                Fn=complex([0,real(Fn),fliplr(real(Fn))],[0,imag(Fn),-fliplr(imag(Fn))]);
                PPn=ifft(Fn);
                PPn=fftshift(PPn);
                P2p=P2p+PPn;
                %Second Order Sub Correction for Potential
                %%% phi=w*t-k.*t*cw;
                Amnp=zeros(1,numel(Apn));
                Amnm=zeros(1,numel(Apn));
                Phii1=zeros(1,numel(Apn));
                Phii2=zeros(1,numel(Apn));
                for j=1:numel(Amnp)
                    if i-j>0
                        Amnp(i-j)=a(i)*a(j)*Am(i,j).*cosh(abs(k(i)-k(j)).*(z+h))./cosh(abs(k(i)-k(j))*h);
                        Phii1(i-j)=phi(i)-phi(j);
                    end
                    if i-j<0
                        Amnm(j-i)=a(i)*a(j)*Am(i,j).*cosh(abs(k(i)-k(j)).*(z+h))./cosh(abs(k(i)-k(j))*h);
                        Phii2(j-i)=-phi(i)+phi(j);
                    end
                end


                Amnp=Amnp*numel(Amnp);
                % Fn=complex(Amn.*cos(Phii),Amn.*sin(Phii));
                Fn=Amnp.*exp(1i*(Phii1-0.5*pi));
                Fn=complex([0,real(Fn),fliplr(real(Fn))],[0,imag(Fn),-fliplr(imag(Fn))]);
                PPn=ifft(Fn);
                PPn=fftshift(PPn);
                P2m=P2m+PPn;

                Amnm=Amnm*numel(Amnm);
                Fn=Amnm.*exp(1i*(Phii2+0.5*pi));
                Fn=complex([0,real(Fn),fliplr(real(Fn))],[0,imag(Fn),-fliplr(imag(Fn))]);
                PPn=ifft(Fn);
                PPn=fftshift(PPn);
                P2m=P2m+PPn;
            end
            P2m=P2m/2;
            P2p=P2p/2;
            %                 PP2m(X_index)=P2m(X_index);
            %                 PP2p(X_index)=P2p(X_index);
            %             end
            %             PP_2=PP_linear+PP2m+PP2p;
            PP_2=PP_linear+P2m+P2p;
            N=length(PP_20);
            PP_20(1:round(0.1*N))=0;
            PP_20(round(0.9*N):end)=0;
            PP_2=PP_2+PP_20;
            %Padding for OW3D
            XX_linear=asymmetric_padding(L,XX_linear);
            PP_linear=asymmetric_padding(L,PP_linear);
            XX_2=asymmetric_padding(L,XX_2);
            PP_2=asymmetric_padding(L,PP_2);

            if numel(XX_linear)~=L
                XX_linear=XX_linear(1:end-1);
                XX_2=XX_2(1:end-1);
                PP_linear=PP_linear(1:end-1);
                PP_2=PP_2(1:end-1);
            end


            %             end
            %%%%%%%%%%Higher order correction
            %Third order
            heta=-imag(hilbert(XX_linear));
            D31=(XX_linear.^2+heta.^2).*XX_linear;
            D33=(XX_linear.^2-3*heta.^2).*XX_linear;
            C=sech(2*kp*h);
            S31=-(3*(1+3*C+3*C^2+2*C^3))/(8*(1-C)^3);
            S33=-S31;
            kp_33=3*kp;
            eta_31=(S31*D31)*kp^2;
            eta_33=(S33*D33)*kp^2;
            XX_33=XX_2+eta_31+eta_33;
            %Third order potential
            pot_31=g/(wp)*-imag(hilbert(eta_31));

            w_p33=sqrt(kp_33*g*tanh(kp_33*h));
            pot_33=g/(w_p33)*-imag(hilbert(eta_33));
            PP_33=PP_2+pot_33+pot_31;
            PP_33=asymmetric_padding(L,PP_33);
            XX_33=asymmetric_padding(L,XX_33);
            %Fourth order
            heta=-imag(hilbert(XX_linear));
            D42=(XX_linear.^2+heta.^2).*(XX_linear.^2-heta.^2);
            D44=(XX_linear.^2-heta.^2).^2-(2.*XX_linear.*heta).^2;
            C=sech(2*kp*h);
            S42=coth(kp*h)*(6-26*C-182*C^2-204*C^3-25*C^4+26*C^5)/((6*(3+2*C)*(1-C)^4));
            S44=coth(kp*h)*(24+92*C+122*C^2+66*C^3+67*C^4+34*C^5)/((24*(3+2*C)*(1-C)^4));
            kp_44=4*kp;
            eta_42=(S42*D42)*kp^3;
            eta_44=(S44*D44)*kp^3;
            XX_44=XX_33+eta_42+eta_44;
            %Third order potential
            pot_42=g/(wp)*-imag(hilbert(eta_42));

            w_p44=sqrt(kp_44*g*tanh(kp_44*h));
            pot_44=g/(w_p44)*-imag(hilbert(eta_44));
            PP_44=PP_33+pot_42+pot_44;
            PP_44=asymmetric_padding(L,PP_44);
            XX_44=asymmetric_padding(L,XX_44);
            %Fifth order
            heta=-imag(hilbert(XX_linear));
            D51=(XX_linear.^2+heta.^2).^2.*XX_linear;
            D53=(XX_linear.^2+heta.^2).*XX_linear.*(XX_linear.^2-3*heta.^2);
            D55=XX_linear.*((XX_linear.^2-heta.^2).^2-(2.*XX_linear.*heta).^2)-4.*heta.^2.*XX_linear.*(XX_linear.^2-heta.^2);
            C=sech(2*kp*h);
            S53 = (9*(132 + 17*C - 2216*C^2 - 5897*C^3 - 6292*C^4 - 2687*C^5 + 194*C^6 + 467*C^7 + 82*C^8)) / (128*(3 + 2*C)*(4 + C)*(1 - C)^6);
            S55 = (5*(300 + 1579*C + 3176*C^2 + 2949*C^3 + 1188*C^4 + 675*C^5 + 1326*C^6 + 827*C^7 + 130*C^8)) / (384*(3 + 2*C)*(4 + C)*(1 - C)^6);
            S51=-(S53+S55);
            kp_55=5*kp;
            eta_53=(S53*D53)*kp^4;
            eta_55=(S55*D55)*kp^4;
            eta_51=(S51*D51)*kp^4;
            XX_55=XX_44+eta_53+eta_55+eta_51;
            %Third order potential
            pot_53=g/(wp)*-imag(hilbert(eta_53));

            w_p55=sqrt(kp_55*g*tanh(kp_55*h));
            pot_55=g/(w_p55)*-imag(hilbert(eta_55));
            pot_51=g/(w_p55)*-imag(hilbert(eta_51));
            PP_55=PP_44+pot_53+pot_55+pot_51;
            PP_55=asymmetric_padding(L,PP_55);
            XX_55=asymmetric_padding(L,XX_55);



            % Visulize
            xx=1:numel(XX_linear);
            xx=xx*dx;
            xx=xx-xx(end)/2;
            switch ORDER
                case 1
                    figure(etah);
                    hold on
                    plot(xx,XX_linear);
                    figure(poth)
                    hold on
                    plot(xx,PP_linear);
                case 2

                    figure(etah);
                    hold on
                    plot(xx,XX_2);
                    figure(poth)
                    hold on
                    plot(xx,PP_2);
                case 3
                    figure(etah);
                    hold on
                    plot(xx,XX_33);
                    figure(poth)
                    hold on
                    plot(xx,PP_33);
            end

            grid on;
            box on;
            grid minor;

            %             for ORDER=1:5
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Write the outputs files
            %Write init
            write_path=fullfile(pwd,"test2",sprintf("Order_%d_Alpha_%.1f_Akp_%.3d_phi_%d",ORDER,Alpha,round(Akp*100),phi_shift));
            if exist(write_path, 'dir')
                disp('目录已存在。');
            else
                disp('目录不存在，正在创建目录。');
                % 创建目录
                mkdir(write_path);
                disp(['目录已创建：', write_path]);
            end
            file_name=fullfile(write_path,'OceanWave3D.init');

            switch ORDER
                case 1
                    eta_S22=XX_linear;
                    pot_S22=PP_linear;
                case 2
                    eta_S22=XX_2;
                    pot_S22=PP_2;
                case 3
                    eta_S22=XX_33;
                    pot_S22=PP_33;
                case 4
                    eta_S22=XX_44;
                    pot_S22=PP_44;
                case 5
                    eta_S22=XX_55;
                    pot_S22=PP_55;

            end
            nx=length(xx);

            ny=1;

            dx=xx(2)-xx(1);

            dy=10;


            Lx=max(xx)-min(xx);

            Ly=ny*dy;

            % dt=0.1;

            fileid=fopen(file_name,'w');

            fprintf(fileid,' H=%f nx=%d ny=%d dx=%f dy=%f akp=%f shift=%f' ,max(eta_S22),nx,ny,dx,dy,Akp,phi_shift/pi);

            fprintf(fileid,'\n%12e %12e %d %d %12e',Lx,Ly,nx,ny,dt);

            for ry=1:ny
                for rx=1:nx
                    fprintf(fileid,'\n%12e %12e',eta_S22(rx),pot_S22(rx));
                end
            end

            fclose(fileid);

            %Write inp
            file_name=fullfile(write_path,'OceanWave3D.inp');
            fileid=fopen(file_name,'w');

            fprintf(fileid,'Semi-Gaussian Modified in %s <-\n',datestr(now,0));
            % fprintf(fileid,'-1 2 %1.2f <-\n',Relative_acceleration)
            if RA_on==1
                fprintf(fileid,'-1 2 0.4 <-\n');
            else
                fprintf(fileid,'-1 2 <-\n');
            end
            fprintf(fileid,'%d 0. %d %d 1 25 0 0 1 1 1 1 <-\n',L*dx,h,L);
            fprintf(fileid,'3 0 3 1 1 1 <-\n');
            % fprintf(fileid,'3 3 3 1 1 1 <-\n')
            fprintf(fileid,'%d %f 1 0. 1 <-\n',round(N_steps)+1,dt);
            fprintf(fileid,'9.81 <-\n');
            fprintf(fileid,'1 3 0 55 1e-6 1e-6 1 V 1 1 20 <-\n');
            fprintf(fileid,'0.05 1.00 1.84 2 0 0 1 6 32 <-\n');
            fprintf(fileid,'8 1 <-\n');
            %     fprintf(fileid,'4 20 1 1 <-\n')
            %     fprintf(fileid,'2 4096 1 1 1 1 1 %d 1 <-\n',round(N_steps))

            fprintf(fileid,'1 0 <-\n');
            fprintf(fileid,'0 6 10 0.08 0.08 0.4 <-\n');
            fprintf(fileid,'0 8. 3 X 0.0 <-\n');
            fprintf(fileid,'0 0 <-\n');
            fprintf(fileid,'0 2.0 2 0 0 1 0 <-\n');
            fprintf(fileid,'0 <-\n');
            fprintf(fileid,'33  8. 2. 80. 20. -1 -11 100. 50. run06.el 22.5 1.0 3.3 <-\n');
            fclose(fileid);


            file_name=fullfile(write_path,'OW_readme.txt');
            fileid=fopen(file_name,'w');
            fprintf(fileid,'Check phase separation T=%1.2fTp H=%f dx=%f, dt=%f\n',t(1)./Tp,max(eta_S22),dx,dt);
            fprintf(fileid,'A=%f,Akp=%f, kp=%f, lambda=%f, Tp=%f\n',A,Akp,kp,lambda,Tp);
            fprintf(fileid,'lambda=%f dx, T=%f dt, cw=%f, CFL=%f, kd=%f\n',lambda/dx,Tp/dt,cw,CFL,kp*h);
            fprintf(fileid,'Last modified in %s',datestr(now,0));
            fclose(fileid);
            fprintf(sprintf('Akp=%.2f,Alpha=%.1f\n',Akp,Alpha));
        end
    end
end
% end
% end
function padded_array = asymmetric_padding(L,array)
%This code is to perform asymmertic zero padding to the original array
% 计算总共需要填充的元素数量
totalPadding = L - numel(array);

% 根据比例计算左右两侧的填充数量
leftRatio = 1;
rightRatio = 20;
totalRatio = leftRatio + rightRatio;

% 计算左侧和右侧应填充的元素数
leftPad = round(totalPadding * (leftRatio / totalRatio));
rightPad = totalPadding - leftPad;

% 对 XX_linear 进行非对称填充
padded_array = padarray(array, [0 leftPad], 0, 'pre');
padded_array = padarray(padded_array, [0 rightPad], 0, 'post');

for i=leftPad:-1:2
    padded_array(i-1)=padded_array(i)/2;
end

for i=1:rightPad
    padded_array(leftPad+length(array)+i)=padded_array(leftPad+length(array)+i-1)/2;
end

end





function resultArray=a20array(originalArray)
% 获取原始数组的长度
originalLength = length(originalArray);

% 创建一个新的数组，长度是原始数组的两倍
newArray = zeros(1, 2 * originalLength);

% 在新数组中插入原始数组元素和0元素
newArray(1:2:end) = originalArray;

% 将新数组的长度减半
resultArray = newArray(1:originalLength);

end

