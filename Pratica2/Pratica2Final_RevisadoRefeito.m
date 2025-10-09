%Aluno: Marcos Mauricio da Cunha de Morais

clear
clc
close all;
pkg load signal

try
  disp("Carregando audio. ");
  pkg load audio
catch
  disp("Audio nao encontrado. Instalando audio. ");
  pkg install -forge audio
  pkg load audio
end

%variaveis globais
numa=[1 0 0.49 0 0 0 0.2401 0 -0.0576 0 -0.282 0 -0.0138]
dena=[1]
av=[0.7 0.9]
Lv=[1 4 10]

disp("Qual Exercicio quer simular? Digite o Numero de 1 a 7: ");
option=input("")

switch option
  case 1
    [h,w]=freqz(numa,dena);

    figure(1);
    subplot(2,1,1);
    plot(w/pi, 20*log10(abs(h)));
    title("Magnitude");
    xlabel("Frequencia Normalizada");
    ylabel("Magnitude (dB)");
    grid on;

    subplot(2,1,2);
    plot(w/pi, unwrap(angle(h))*180/pi);
    title("Resposta em Frequencia: Fase");
    xlabel("Frequencia Normalizada");
    ylabel("Fase (graus)");
    grid on;

    figure(2);
    zeroshorig=roots(numa);
    poloshorig=zeros(length(zeroshorig),1);
    zplane(zeroshorig,poloshorig);
    title("Diagrama de Polos e Zeros");

  case 2
    load handel.mat
    y_filtrado=filter(numa,dena,y);

    spectrumAnalyzer(y_filtrado,Fs);
    title("Espectro da Resposta de H(z) ao audio handel.mat");
    sound(y_filtrado,Fs);

  case 3
    invnuma=dena
    invdena=numa

    [hinv,w]=freqz(invnuma,invdena);

    load handel.mat

    %NOTA: depois de 2 cascos azuis e cair na rainbow road pra 12 lugar, to aqui dnv

    y_filtrado=filter(numa,dena,y);
    y_recuperado=filter(invnuma,invdena,y_filtrado);

    figure(1);
    subplot(2,1,1);
    plot(w/pi, 20*log10(abs(hinv)));
    title("Resposta em Frequencia Inversa");
    xlabel("Frequencia Normalizada");
    ylabel("Magnitude (dB)");
    grid on;

    subplot(2,1,2);
    plot(w/pi, unwrap(angle(hinv))*180/pi);
    xlabel("Frequencia Normalizada");
    ylabel("Fase (graus)");
    grid on;

    figure(2);
    polosinv=roots(invdena);
    zerosinv=zeros(length(polosinv),1);
    zplane(zerosinv,polosinv);
    title("Diagrama de Polos e Zeros Inversos");

    nmin=min(length(y),length(y_recuperado));
    mse=mean((y(1:nmin)-y_recuperado(1:nmin)).^2);

    disp("Erro quadratico medio entre original e recuperado: ");
    disp(mse);

    spectrumAnalyzer(y,Fs);
    title("Espectro do Sinal Original");

    spectrumAnalyzer(y_recuperado,Fs);
    title("Espectro do Sinal Recuperado");

    sound(y_recuperado,Fs);

  case 4
    figc=1;
    plotind=1;

    for a=av
      for L=Lv
        numera=[1 zeros(1,L-1) -1];
        denomi=[1 zeros(1,L-1) -a];

        [h,w]=freqz(numera,denomi);

        figure(figc);
        figc=figc+1;

        subplot(2, 1, 1);
        plot(w/pi, 20*log10(abs(h)));
        title(sprintf("Magnitude (a=%.1f, L=%d)", a, L));
        xlabel('Frequencia Normalizada');
        ylabel("Magnitude (dB)");
        grid on;

        subplot(2, 1, 2);
        plot(w/pi, unwrap(angle(h))*180/pi);
        title("Fase");
        xlabel("Frequencia Normalizada");
        ylabel("Fase (graus)");
        grid on;

        zeroscomb=roots(numera);
        poloscomb=roots(denomi);

        figure(figc);
        figc=figc+1;

        zplane(zeroscomb, poloscomb);
        title(sprintf("Polos e Zeros (a=%.1f, L=%d)", a, L));

        plotind=plotind+1;
      endfor
    endfor

  case 5
    load handel.mat
    for a=av
      for L=Lv
        numera=[1 zeros(1,L-1) -1];
        denomi=[1 zeros(1,L-1) -a];
        yout=filter(numera,denomi,y);

        spectrumAnalyzer(yout,Fs);
        title(sprintf("Espectro da Resposta (a=%.1f, L=%d)", a, L));
        sound(yout,Fs);
      endfor
    endfor

  case 6
    load handel.mat
    figc=0;

    for a=av
      for L=Lv
        numera=[1 zeros(1,L-1) -1];
        denomi=[1 zeros(1,L-1) -a];
        numinv=denomi
        deninv=numera

        %NOTA: que vontade de jogar super mario galaxy 2

        y_filtrado=filter(numera,denomi,y);
        y_recuperado=filter(numinv,deninv,y_filtrado);

        [hinv,w]=freqz(numinv,deninv);

        figure(figc);
        figc=figc+1;

        subplot(2, 1, 1);
        plot(w/pi, 20*log10(abs(hinv)));
        title(sprintf("Magnitude Inversa (a=%.1f, L=%d)", a, L));
        xlabel("Frequencia Normalizada");
        ylabel("Magnitude (dB)");
        grid on;

        polosinv=roots(deninv);
        zerosinv=roots(numinv);
        figure(figc);
        figc=figc+1;

        zplane(zerosinv,polosinv);
        title(sprintf("Polos e Zeros Inversos (a=%.1f, L=%d)", a, L));

        nmin=min(length(y),length(y_recuperado));
        mse=mean((y(1:nmin)-y_recuperado(1:nmin)).^2);
        disp("Erro quadratico medio entre original e recuperado: ");
        disp(mse);

        %sound(y_recuperado,Fs);
      endfor
    endfor

  case 7
    norder=[3 5 7]
    hlength=50

    figc=1;
    numinvq3=dena;
    deninvq3=numa;

    %NOTA: terminar amanha

    [hq3,dummy]=deconv([numinvq3,zeros(1,hlength)],deninvq3);
    [href,wref]=freqz([numinvq3, zeros(1,hlength)],deninvq3);

    for N=norder
      num_fir=hq3(1:N+1);
      den_fir=[1]
      [h_fir,w_fir]=freqz(num_fir,den_fir);

      figure(figc);
      figc=figc+1;

      plot(wref/pi, 20*log10(abs(href)), "b", "LineWidth", 2);
      hold on;
      plot(w_fir/pi, 20*log10(abs(h_fir)), "r--", "LineWidth", 1.5);
      hold off;

      title(sprintf("Q3: FIR Ordem %d vs IIR - Magnitude", N));
      xlabel("Frequencia Normalizada");
      ylabel("Magnitude (dB)");
      legend("IIR (Referencia)", sprintf("FIR Ordem %d", N));
      grid on;
    endfor

    for a=av
      for L=Lv
        numera=[1 zeros(1,L-1) -1];
        denomi=[1 zeros(1,L-1) -a];
        numinvq6=denomi;
        deninvq6=numera;

        num_deconv=[numinvq6,zeros(1,hlength+length(deninvq6) - 1)];
        [hq6,dummy]=deconv(num_deconv,deninvq6);
        [h_ref,w_ref]=freqz(numinvq6,deninvq6);

        for N=norder
          num_fir=hq6(1:N+1);

          [h_fir,w_fir]=freqz(num_fir,[1]);

          figure(figc);
          figc=figc+1;

          plot(w_ref/pi, 20*log10(abs(h_ref)), "b", "LineWidth", 2);
          hold on;
          plot(w_fir/pi, 20*log10(abs(h_fir)), "r--", "LineWidth", 1.5);
          hold off;

          title(sprintf("Q6: Comb a=%.1f, L=%d - FIR Ordem %d vs IIR", a, L, N));
          xlabel("Frequencia Normalizada");
          ylabel("Magnitude (dB)");
          legend("IIR (Referencia)", sprintf("FIR Ordem %d", N));
          grid on;

          %NOTA: 21 graficos e brincadeira, pra quebrar a banca do silvio santos

        endfor
      endfor
    endfor

  otherwise
    disp("Opcao invalida. Terminando o Programa.");
    pause(2);
end
