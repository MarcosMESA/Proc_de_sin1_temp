%Aluno: Marcos Mauricio da Cunha de Morais
%Aula Lista 3

clc;
close all;
clear;

try
  disp("Carregando audio. ");
  pkg load audio
catch
  disp("Audio nao encontrado. Instalando audio. ");
  pkg install -forge audio
  pkg load audio
end

try
  disp("Carregando signal. ");
  pkg load signal
catch
  disp("Signal nao encontrado. Instalando signal. ");
  pkg install -forge signal
  pkg load signal
end

N=1024;

disp("Qual Exercicio quer simular? Digite o Numero de 1 a 5: ");
option=input("");

switch option
  case 1
    wc_lp=0.25; %corte passa baixa
    wc_hp=0.35; %do passa alta
    T=1;

    %passa baixa
    Wc_lp=2/T*tan(wc_lp*pi/2);
    [b_lp,a_lp]=butter(2,wc_lp,'low');

    %passa alta
    Wc_hp=2/T*tan(wc_hp*pi/2);
    [b_hp,a_hp]=butter(2,wc_hp,'high');

    [H_lp,w]=freqz(b_lp,a_lp,N);
    [H_hp,w]=freqz(b_hp,a_hp,N);

    H_combined=H_lp+H_hp;

    figure(1);
    subplot(2,1,1);

    y_limits=[-40 5];
    fill([0.25 0.35 0.35 0.25],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],[0.9 0.9 0.8],'FaceAlpha',0.3,'EdgeColor','none');
    hold on;

    plot(w/pi,20*log10(abs(H_lp)),'b','LineWidth',1.5);
    plot(w/pi,20*log10(abs(H_hp)),'r','LineWidth',1.5);
    plot(w/pi,20*log10(abs(H_combined)),'g','LineWidth',2);

    grid on;
    xlabel('Frequencia Normalizada (×π rad/amostra)');
    ylabel('Magnitude (dB)');
    title('Resposta dos Filtros Individuais e Combinado');
    legend('Passa-Baixas (ω_c=0.25)','Passa-Altas (ω_c=0.35)','Rejeita-Faixa','Location','northeast');
    xlim([0 1]);
    ylim(y_limits);

    subplot(2,1,2);
    y_limits=[0 1.2];
    fill([0.25 0.35 0.35 0.25],[y_limits(1) y_limits(1) y_limits(2) y_limits(2)],[0.9 0.9 0.8],'FaceAlpha',0.3,'EdgeColor','none');
    hold on;

    plot(w/pi,abs(H_combined),'k','LineWidth',2);

    grid on;
    xlabel('Frequência Normalizada (×π rad/amostra)');
    ylabel('Magnitude (linear)');
    title('Resposta do Filtro Rejeita-Faixa (linear)');
    xlim([0 1]);
    ylim(y_limits);

    fprintf('\nPassa-Baixas (2 ordem):\n');
    fprintf('b_lp=[%f,%f,%f]\n',b_lp);
    fprintf('a_lp=[%f,%f,%f]\n',a_lp);
    fprintf('\nPassa-Altas (2 ordem):\n');
    fprintf('b_hp=[%f,%f,%f]\n',b_hp);
    fprintf('a_hp=[%f,%f,%f]\n',a_hp);

    %NOTA: tentei usar uistack(), mas so tem pra matlab, pq vi num exemplo na net

  case 2
    f_notch=10;
    ws=200;

    fs_hz=ws/(2*pi);
    f_nyquist=fs_hz/2;
    f_norm=f_notch/f_nyquist;
    w_notch=f_norm*pi;
    r=0.9;
    b_notch=[1,-2*cos(w_notch),1];
    a_notch=[1,-2*r*cos(w_notch),r^2];
    b_notch=b_notch/sum(b_notch);

    [H,w]=freqz(b_notch,a_notch,N);

    figure(1);

    subplot(2,1,1);
    plot(w/pi,20*log10(abs(H)),'b','LineWidth',2);
    grid on;
    xlabel('Frequencia Normalizada (×π rad/amostra)');
    ylabel('Magnitude (dB)');
    title(sprintf('Filtro Notch de 2 Ordem - Eliminacao de %d Hz',f_notch));
    xlim([0 1]);

    hold on;
    plot([w_notch/pi,w_notch/pi],ylim,'r--','LineWidth',1);
    text(w_notch/pi+0.02,-20,sprintf('10 Hz\n(ω=%.3fπ)',w_notch/pi),'Color','r','FontWeight','bold');

    subplot(2,1,2);
    plot(w/pi,abs(H),'b','LineWidth',2);
    grid on;
    xlabel('Frequencia Normalizada (×π rad/amostra)');
    ylabel('Magnitude (linear)');
    title('Resposta em Magnitude (linear)');
    xlim([0 1]);
    ylim([0 1.1]);

    hold on;
    plot([w_notch/pi,w_notch/pi],ylim,'r--','LineWidth',1);

    fprintf('Coeficientes:\n');
    fprintf('b=[%.6f,%.6f,%.6f]\n',b_notch);
    fprintf('a=[%.6f,%.6f,%.6f]\n',a_notch);

  case 3
    N=4;
    M=N;
    w1=pi/6;
    w2=pi/3;
    h_ideal=zeros(1,M+1);

    for n=0:M
      m=n-M/2;
      if m==0
        h1=w1/pi;
        h2=w2/pi;
      else
        h1=sin(w1*m)/(pi*m);
        h2=sin(w2*m)/(pi*m);
      endif

      h_ideal(n+1)=h1+h2;
    endfor

    w_hamming=zeros(1,M+1);

    for n=0:M
      w_hamming(n+1)=0.54-0.46*cos(2*pi*n/M);
    endfor

    h_fir=h_ideal.*w_hamming;
    N_points=512;
    [H,w]=freqz(h_fir,1,N_points);

    H_ideal=zeros(size(w));
    for i=1:length(w)
      if w(i)<w1
        H_ideal(i)=2;
      elseif w(i)<w2
        H_ideal(i)=1;
      else
        H_ideal(i)=0;
      endif
    endfor

    figure(1);

    subplot(3,1,1);
    plot(w/pi,abs(H),'b-','LineWidth',2);
    hold on;
    plot(w/pi,H_ideal,'r--','LineWidth',1.5);
    grid on;
    xlabel('Frequencia Normalizada (×π rad/amostra)');
    ylabel('Magnitude');
    title('Resposta em Frequencia - FIR 4 ordem com Janela de Hamming');
    legend('Filtro Prático','Filtro Ideal','Location','northeast');
    xlim([0 1]);
    ylim([-0.2 2.5]);

    plot([w1/pi,w1/pi],[-0.2 2.5],'k:','LineWidth',1);
    plot([w2/pi,w2/pi],[-0.2 2.5],'k:','LineWidth',1);
    text(w1/pi/2,2.2,'Ganho=2','HorizontalAlignment','center');
    text((w1/pi+w2/pi)/2,1.5,'Ganho=1','HorizontalAlignment','center');
    text((w2/pi+1)/2,0.5,'Ganho=0','HorizontalAlignment','center');

    subplot(3,1,2);
    stem(0:M,h_fir,'b','LineWidth',1.5,'MarkerFaceColor','b');
    hold on;
    stem(0:M,h_ideal,'r--','LineWidth',1,'MarkerFaceColor','r');
    grid on;
    xlabel('Amostra (n)');
    ylabel('Amplitude');
    title('Coeficientes do Filtro FIR (Resposta ao Impulso)');
    legend('FIR Pratico (Hamming)','FIR Ideal','Location','northeast');

    subplot(3,1,3);
    stem(0:M,w_hamming,'g','LineWidth',1.5,'MarkerFaceColor','g');
    grid on;
    xlabel('Amostra (n)');
    ylabel('Amplitude');
    title('Janela de Hamming');
    ylim([0 1.1]);

    gain_dc=sum(h_fir);
    fprintf('Ganho em DC (ω=0): H(0) = sum(h[n]) = %.6f (ideal: 2.0000)\n',gain_dc);

    idx_w1=round(w1*N_points/pi);
    gain_band1=mean(abs(H(1:idx_w1)));
    fprintf('Ganho medio em [0, π/6): %.6f (ideal: 2.0000)\n',gain_band1);

    idx_w2=round(w2*N_points/pi);
    gain_band2=mean(abs(H(idx_w1+1:idx_w2)));
    fprintf('Ganho medio em [π/6, π/3): %.6f (ideal: 1.0000)\n',gain_band2);

    gain_band3=mean(abs(H(idx_w2+1:end)));
    fprintf('Ganho medio em [π/3, π]: %.6f (ideal: 0.0000)\n',gain_band3);

    impulse_input=[1,zeros(1,20)];
    impulse_response=filter(h_fir,1,impulse_input);

    for i=1:min(10,length(impulse_response))
      fprintf('h[%d] = %.6f\n',i-1,impulse_response(i));
    end

    figure(2);
    subplot(2,1,1);
    plot(w/pi,angle(H),'b-','LineWidth',1.5);
    grid on;
    xlabel('Frequencia Normalizada (×π rad/amostra)');
    ylabel('Fase (radianos)');
    title('Resposta de Fase do Filtro FIR');
    xlim([0 1]);

    subplot(2,1,2);
    plot(w/pi,unwrap(angle(H)),'r-','LineWidth',1.5);
    grid on;
    xlabel('Frequencia Normalizada (×π rad/amostra)');
    ylabel('Fase (radianos - desenrolada)');
    title('Fase Desenrolada');
    xlim([0 1]);

  case 4
    N=20;
    M=N;
    ws=2*pi;
    wc=pi/4;

    %a)
    h_ideal=zeros(1,M+1);

    for n=0:M
      m=n-M/2;
      if m==0
        h_ideal(n+1)=wc/pi;
      else
        h_ret=sin(wc*m)/(pi*m);
        h_ideal(n+1)=h_ret^2;
      endif
    endfor

    h_ideal=h_ideal / sum(h_ideal);

    %b)
    w_triangular=zeros(1,M+1);

    for n=0:M
      if n<=M/2
        w_triangular(n+1)=2*n/M;
      else
        w_triangular(n+1)=2-2*n/M;
      endif
    endfor

    h_fir=h_ideal.*w_triangular;
    h_fir = h_fir / sum(h_fir);

    N_points=512;
    [H,w]=freqz(h_fir,1,N_points);

    H_ideal=zeros(size(w));
    for i=1:length(w)
      freq=w(i);
      if abs(freq)<=2*wc
        H_ideal(i)=1-abs(freq)/(2*wc);
      else
        H_ideal(i)=0;
      endif
    endfor

    figure(1);

    subplot(3,1,1);
    plot(w/pi,abs(H),'b-','LineWidth',2);
    hold on;
    plot(w/pi,H_ideal,'r--','LineWidth',1.5);
    grid on;
    xlabel('Frequencia Normalizada (×π rad/amostra)');
    ylabel('Magnitude |H(e^{jω})|');
    title('Resposta em Frequencia - FIR 4 ordem (Janela Triangular)');
    legend('Filtro Pratico','Filtro Ideal','Location','northeast');
    xlim([0 1]);
    ylim([-0.1 1.2]);

    plot([2*wc/pi,2*wc/pi],[-0.1 1.2],'k:','LineWidth',1);
    text(2*wc/pi+0.05,0.5,sprintf('2ω_c=π/2≈%.3fπ',2*wc/pi));

    subplot(3,1,2);
    stem(0:M,h_fir,'b','LineWidth',1.5,'MarkerFaceColor','b');
    hold on;
    stem(0:M,h_ideal,'r--','LineWidth',1,'MarkerFaceColor','r');
    grid on;
    xlabel('Amostra (n)');
    ylabel('Amplitude');
    title('Coeficientes do Filtro FIR');
    legend('FIR Pratico','FIR Ideal','Location','northeast');

    subplot(3,1,3);
    stem(0:M,w_triangular,'g','LineWidth',1.5,'MarkerFaceColor','g');
    grid on;
    xlabel('Amostra (n)');
    ylabel('Amplitude');
    title('Janela Triangular (Bartlett)');
    ylim([0 1.1]);


    %c)

    gain_dc=sum(h_fir);
    fprintf('Ganho em DC (ω=0): %.6f (ideal: 1.0000)\n', gain_dc);

    [H,freq]=freqz(h_fir,1,512);
    [~,idx]=min(abs(freq-wc));
    H_wc=H(idx);
    fprintf('Ganho em ωc=π/4: %.6f (ideal: 0.5000)\n', abs(H_wc));

    omega_2wc=2*wc;
    H_2wc=0;
    for n=0:length(h_fir)-1
      H_2wc=H_2wc+h_fir(n+1)*exp(-1j*omega_2wc*n);
    endfor
    fprintf('Ganho em 2ωc=π/2: %.6f (ideal: 0.0000)\n', abs(H_2wc));

    %NOTA: eu acho q isso ta errado, so acho kkkkkkkk tentei consertar e melhorou um pouco, mas nao tanto

  case 5
    N=20;
    M=N;
    w1=pi/4;
    w2=pi/2;

    %a)
    h_ideal=zeros(1,M+1);

    for n=0:M
      m=n-M/2;
      if m==0
        h_ideal(n+1)=(w2-w1)/pi;
      else
        h1=sin(w1*m)/(pi*m);
        h2=sin(w2*m)/(pi*m);
        h_ideal(n+1)=h2-h1;
      endif
    endfor

    w_center=(w1+w2)/2;
    H_center=0;

    for n=0:M
      H_center=H_center+h_ideal(n+1)*exp(-1j*w_center*n);
    endfor

    h_ideal=h_ideal/abs(H_center);

    %b)
    w_hann=zeros(1,M+1);

    for n=0:M
      w_hann(n+1)=0.5-0.5*cos(2*pi*n/M);
    endfor

    h_fir=h_ideal.*w_hann;
    H_center=0;

    for n=0:M
      H_center=H_center+h_fir(n+1)*exp(-1j*w_center*n);
    endfor

    h_fir=h_fir/abs(H_center);

    N_points=512;
    [H,w]=freqz(h_fir,1,N_points);
    H_ideal=zeros(size(w));

    for i=1:length(w)
      freq=w(i);
      if freq<w1
        H_ideal(i)=freq/w1;
      elseif freq<w2
        H_ideal(i)=1;
      else
        H_ideal(i)=0;
      endif
    endfor

    figure(1);

    subplot(3,1,1);
    plot(w/pi,abs(H),'b-','LineWidth',2);
    hold on;
    plot(w/pi,H_ideal,'r--','LineWidth',1.5);
    grid on;
    xlabel('Frequencia Normalizada (×π rad/amostra)');
    ylabel('Magnitude |H(e^{jω})|');
    title('Resposta em Frequencia - FIR 4 ordem (Janela de Hann)');
    legend('Filtro Pratico','Filtro Ideal','Location','northeast');
    xlim([0 1]);
    ylim([-0.1 1.2]);

    plot([w1/pi,w1/pi],[-0.1 1.2],'k:','LineWidth',1);
    plot([w2/pi,w2/pi],[-0.1 1.2],'k:','LineWidth',1);
    text(w1/(2*pi),0.6,'Rampa','HorizontalAlignment','center');
    text((w1/pi+w2/pi)/2,1.1,'Ganho=1','HorizontalAlignment','center');
    text((w2/pi+1)/2,0.5,'Ganho=0','HorizontalAlignment','center');

    subplot(3,1,2);
    stem(0:M,h_fir,'b','LineWidth',1.5,'MarkerFaceColor','b');
    hold on;
    stem(0:M,h_ideal,'r--','LineWidth',1,'MarkerFaceColor','r');
    grid on;
    xlabel('Amostra (n)');
    ylabel('Amplitude');
    title('Coeficientes do Filtro FIR');
    legend('FIR Pratico','FIR Ideal','Location','northeast');

    subplot(3,1,3);
    stem(0:M,w_hann,'g','LineWidth',1.5,'MarkerFaceColor','g');
    grid on;
    xlabel('Amostra (n)');
    ylabel('Amplitude');
    title('Janela de Hann (Hanning)');
    ylim([0 1.1]);

    %c)
    H_w1=0;
    for n=0:length(h_fir)-1
      H_w1=H_w1+h_fir(n+1)*exp(-1j*w1*n);
    endfor
    fprintf('Ganho em w1=π/4: %.6f (ideal: 1.0000)\n', abs(H_w1));

    H_wcenter=0;
    for n=0:length(h_fir)-1
      H_wcenter=H_wcenter+h_fir(n+1)*exp(-1j*w_center*n);
    endfor
    fprintf('Ganho em (w1+w2)/2: %.6f (ideal: 1.0000)\n', abs(H_wcenter));

    H_w2=0;
    for n=0:length(h_fir)-1
      H_w2=H_w2+h_fir(n+1)*exp(-1j*w2*n);
    endfor
    fprintf('Ganho em w2=π/2: %.6f (ideal: 0.0000)\n', abs(H_w2));

    H_0=sum(h_fir);
    fprintf('Ganho em ω=0 (DC): %.6f (ideal: 0.0000)\n', abs(H_0));

    H_pi=0;
    for n=0:length(h_fir)-1
      H_pi=H_pi+h_fir(n+1)*exp(-1j*pi*n);
    endfor
    fprintf('Ganho em ω=π: %.6f (ideal: 0.0000)\n', abs(H_pi));

  otherwise
    disp("Opcao invalida. Terminando o Programa.");
    pause(2);

endswitch
