%Aluno: Marcos Mauricio
%Processamento de Sinais 1 - Aula Prática 1

%Apresentacao

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

disp("Qual Exercicio quer simular? Digite o Numero de 1 a 7: ");
option = input("");

%Exercicios

switch option
  case 1
    fs=44100;
    t=0:(1/fs):5;
    f=[500 5000 10000];
    x=cos(2*pi* f' *t);
    fst_sign=x(1,:);
    snd_sign=x(2,:);
    trd_sign=x(3,:);

    figure;
    set(gcf,'position',[0 0 800 600]);
    hold on;
    plot(t,fst_sign,'linewidth',2);
    plot(t,snd_sign,'linewidth',2);
    plot(t,trd_sign,'linewidth',2);

    spectrumAnalyzer(fst_sign,fs)
    spectrumAnalyzer(snd_sign,fs)
    spectrumAnalyzer(trd_sign,fs)

    sound(fst_sign,fs)
    pause(3)
    sound(snd_sign,fs)
    pause(3)
    sound(trd_sign,fs)
    pause(3)

  case 2
    fs=44100;
    t=0:(1/fs):5;
    f=[500 10000];

    rsull=chirp(t,f(1),5,f(2),'linear');
    rsulq=chirp(t,f(1),5,f(2),'quadratic');
    rsullog=chirp(t,f(1),5,f(2),'logarithmic');

    spectrumAnalyzer(rsull,fs)
    spectrumAnalyzer(rsulq,fs)
    spectrumAnalyzer(rsullog,fs)

    sound(rsull,fs)
    pause(3)
    sound(rsulq,fs)
    pause(3)
    sound(rsullog,fs)
    pause(3)

  case 3
    load handel.mat;
    figure;
    set(gcf, 'position', [0 0 1000 600]);
    t=(0:length(y)-1)/Fs;
    plot(t, y);
    xlabel('Tempo (s)', 'fontname', 'times new roman', 'fontsize', 14);
    ylabel('Amplitude', 'fontname', 'times new roman', 'fontsize', 14);
    title('Sinal de áudio - Handel (Domínio do Tempo)', 'fontname', 'times new roman', 'fontsize', 16);
    grid on;

    spectrumAnalyzer(y, Fs);
    title('Espectro do sinal de áudio - Handel', 'fontname', 'times new roman', 'fontsize', 16);
    sound(y, Fs);
    pause(length(y)/Fs + 1);  % Aguardar o término da reprodução + 1 segundo
    sound(y, 2*Fs);
    pause(length(y)/(2*Fs) + 1);
    sound(y, 4*Fs);
    pause(length(y)/(4*Fs) + 1);

  case 4
    load handel.mat;
    M_factors=[2, 4, 8];
    for i = 1:length(M_factors)
      M=M_factors(i);
      y_upsampled=upsample(y, M);
      Fs_new=Fs*M;
      figure;
      spectrumAnalyzer(y_upsampled, Fs_new);
      title(sprintf('Espectro após upsample (M=%d)', M));
      sound(y_upsampled,Fs_new);
      pause(length(y_upsampled)/Fs_new + 1);
    endfor

    for i = 1:length(M_factors)
      M=M_factors(i);
      y_resampled=resample(y,M,1);
      Fs_new=Fs*M;
      figure;
      spectrumAnalyzer(y_resampled,Fs_new);
      title(sprintf('Espectro após resample (M=%d)', M));
      sound(y_resampled,Fs_new);
      pause(length(y_resampled)/Fs_new+1);
    endfor

    for i = 1:length(M_factors)
      M=M_factors(i);
      y_downsampled=downsample(y,M);
      Fs_new=Fs/M;
      figure;
      spectrumAnalyzer(y_downsampled,Fs_new);
      title(sprintf('Espectro após downsample (M=%d)', M));
      sound(y_downsampled,Fs_new);
      pause(length(y_downsampled)/Fs_new+1);
    endfor

    for i = 1:length(M_factors)
      M=M_factors(i);
      y_resampled=resample(y,1,M);
      Fs_new=Fs/M;
      figure;
      spectrumAnalyzer(y_resampled,Fs_new);
      title(sprintf('Espectro após resample (1/%d)', M));
      sound(y_resampled,Fs_new);
      pause(length(y_resampled)/Fs_new+1);
    endfor

  case 5
    disp("Nada acontece feijoada.");

  case 6
    audio_files={'h_banheiro.wav','sinal_taca.wav'};
    for i = 1:length(audio_files)
      filename = audio_files{i};
      [y, Fs] = audioread(filename);
      if size(y,2) > 1
        y=mean(y,2);
      endif

      figure;
      set(gcf,'position',[0, 0, 1000, 600]);
      t=(0:length(y)-1)/Fs;
      plot(t,y);
      xlabel('Tempo (s)', 'fontname', 'times new roman', 'fontsize', 14);
      ylabel('Amplitude', 'fontname', 'times new roman', 'fontsize', 14);
      title(sprintf('%s - Domínio do Tempo', filename), 'fontname', 'times new roman', 'fontsize', 16);
      grid on;

      spectrumAnalyzer(y,Fs);
      title(sprintf('Espectro de %s', filename));

      sound(y,Fs);
      pause(length(y)/Fs+1);
    endfor

  case 7
    [y_banheiro, Fs_banheiro]=audioread('h_banheiro.wav');
    [y_taca, Fs_taca]=audioread('sinal_taca.wav');
    if size(y_banheiro, 2) > 1
      y_banheiro = mean(y_banheiro, 2);
    endif
    if size(y_taca, 2) > 1
      y_taca = mean(y_taca, 2);
    endif
    if Fs_banheiro ~= Fs_taca
      y_taca = resample(y_taca, Fs_banheiro, Fs_taca);
      Fs=Fs_banheiro;
    else
      Fs=Fs_banheiro;
    endif

    y_banheiro=y_banheiro/max(abs(y_banheiro));
    y_taca=y_taca/max(abs(y_taca));

    conv_banheiro_taca=conv(y_banheiro, y_taca);
    conv_taca_banheiro=conv(y_taca, y_banheiro);

    figure;
    subplot(3, 1, 1);
    t_banheiro=(0:length(y_banheiro)-1)/Fs;
    plot(t_banheiro, y_banheiro);
    title('Resposta ao Impulso do Banheiro (h_{banheiro}[n])');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    grid on;

    subplot(3, 1, 2);
    t_taca=(0:length(y_taca)-1)/Fs;
    plot(t_taca, y_taca);
    title('Sinal da Taça');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    grid on;

    subplot(3, 1, 3);
    t_conv=(0:length(conv_banheiro_taca)-1)/Fs;
    plot(t_conv, conv_banheiro_taca);
    title('Convolução: h_{banheiro}[n] * sinal_{taça}[n]');
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    grid on;

    figure;
    subplot(3, 1, 1);
    spectrumAnalyzer(y_banheiro, Fs);
    title('Espectro da Resposta ao Impulso do Banheiro');
    subplot(3, 1, 2);
    spectrumAnalyzer(y_taca, Fs);
    title('Espectro do Sinal da Taça');
    subplot(3, 1, 3);
    spectrumAnalyzer(conv_banheiro_taca, Fs);
    title('Espectro da Convolução');

    sound(y_banheiro, Fs);
    pause(length(y_banheiro)/Fs + 1);
    sound(y_taca, Fs);
    pause(length(y_taca)/Fs + 1);
    sound(conv_banheiro_taca, Fs);
    pause(length(conv_banheiro_taca)/Fs + 1);
    sound(conv_taca_banheiro, Fs);
    pause(length(conv_taca_banheiro)/Fs + 1);

  otherwise
    disp("Opcao invalida. Terminando o Programa.");
    pause(2);
end


