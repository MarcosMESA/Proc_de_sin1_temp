%Aluno: Marcos Mauricio da Cunha de Morais
%Aula Pratica 4

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

%variaveis globais
fs=20000;
Ts=1/fs;
r_values=[0.9, 0.95, 0.98];

disp("Qual Exercicio quer simular? Digite o Numero de 1 a 5: ");
option=input("")

switch option
  case 1
    %a) passa baixa
    %NOTA: eu ia usar nome em portugues pras variaveis, mas fica mt "feio"
    fc_low=1000;
    wc_low=2*pi*fc_low/fs;

    colors=['r', 'g', 'b'];

    figure('Name', 'Filtro Passa-Baixas', 'NumberTitle', 'off');
    hold on;

    for idx=1:length(r_values)
      r=r_values(idx);
      b_low=1 - r; %num
      a_low=[1, -r]; %den

      %resposta em freq
      %NOTA: ja foi usado freqz em outra aula, eu acho, entao to usando aqui
      %NoTA: Foi na Aula 2
      [H_low, f_low]=freqz(b_low, a_low, 1024, fs);
      plot(f_low, 20*log10(abs(H_low)), colors(idx), 'LineWidth', 1.5,'DisplayName', ['r=', num2str(r)]);
    endfor

    hold off;
    grid on;
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Filtro Passa-Baixas - fc=1000 Hz');
    legend('show');
    xlim([0, fs/2]);

    %b) passa alta
    fc_high=2000;
    wc_high=2*pi*fc_high/fs;

    figure('Name', 'Filtro Passa-Altas', 'NumberTitle', 'off');
    hold on;

    for idx=1:length(r_values)
      r=r_values(idx);
      b_high=[1, -1] * (1+r)/2;
      a_high=[1, -r];
      [H_high, f_high]=freqz(b_high, a_high, 1024, fs);
      plot(f_high, 20*log10(abs(H_high)), colors(idx), 'LineWidth', 1.5,'DisplayName', ['r=', num2str(r)]);
    endfor

    hold off;
    grid on;
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Filtro Passa-Altas - fc=2000 Hz');
    legend('show');
    xlim([0, fs/2]);

    %c) filtro notch
    %NOTa: Não o do minecraft

    fc_notch=3000;
    wc_notch=2*pi*fc_notch/fs;

    figure('Name', 'Filtro Notch', 'NumberTitle', 'off');
    hold on;

    for idx=1:length(r_values)
      r=r_values(idx);

      %filtro notch de 2 ordem
      cos_w0=cos(wc_notch);
      b_notch=[1, -2*cos_w0, 1];
      a_notch=[1, -2*r*cos_w0, r^2];

      [H_notch, f_notch]=freqz(b_notch, a_notch, 1024, fs);
      plot(f_notch, 20*log10(abs(H_notch)), colors(idx), 'LineWidth', 1.5,'DisplayName', ['r=', num2str(r)]);
    endfor

    hold off;
    grid on;
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Filtro Notch - fc=3000 Hz');
    legend('show');
    xlim([0, fs/2]);
    ylim([-50, 5]);

  case 2
    fc1_bp=500; %corte inferior
    fc2_bp=2000; %corte superior
    wc1_bp=2*pi*fc1_bp/fs;
    wc2_bp=2*pi*fc2_bp/fs;

    colors=['r', 'g', 'b']; %red green blue (nao e pokemon)
    figure('Name', 'Filtro Passa-Faixa', 'NumberTitle', 'off');
    hold on;

    for idx = 1:length(r_values)
      r = r_values(idx);
      cos_wc1=cos(wc1_bp);
      b_bp1=[1, 0, -1];
      a_bp1=[1, -2*r*cos_wc1, r^2];

      cos_wc2= cos(wc2_bp);
      b_bp2=[1, 2, 1]; %ajustar ganho
      a_bp2=[1, -2*r*cos_wc2, r^2];

      %combinar 2 filtros
      b_bp_total=conv(b_bp1, b_bp2);
      a_bp_total=conv(a_bp1, a_bp2);

      [H_bp, f_bp] = freqz(b_bp_total, a_bp_total, 1024, fs);
      H_bp_max = max(abs(H_bp));
      b_bp_total = b_bp_total / H_bp_max;

      [H_bp, f_bp] = freqz(b_bp_total, a_bp_total, 1024, fs); %frequencia normalizada

      plot(f_bp, 20*log10(abs(H_bp)), colors(idx), 'LineWidth', 1.5,'DisplayName', ['r = ', num2str(r)]);
    endfor

    hold off;
    grid on;
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title(['Filtro Passa-Faixa: fc1 = ', num2str(fc1_bp), ' Hz, fc2 = ', num2str(fc2_bp), ' Hz']);
    legend('show');
    xlim([0, fs/2]);
    ylim([-40, 5]);

    hold on; %linhas verticais na freq de corte
    plot([fc1_bp, fc1_bp], ylim, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot([fc2_bp, fc2_bp], ylim, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    hold off;

  case 3
    fc1_bs=1000;
    fc2_bs=2500;
    wc1_bs=2*pi*fc1_bs/fs;
    wc2_bs=2*pi*fc2_bs/fs;

    colors = ['r', 'g', 'b'];
    figure('Name', 'Filtro Rejeita-Faixa', 'NumberTitle', 'off');
    hold on;

    for idx = 1:length(r_values)
      r=r_values(idx);

      cos_wc1=cos(wc1_bs);
      b_notch1=[1, -2*cos_wc1, 1];
      a_notch1=[1, -2*r*cos_wc1, r^2];
      cos_wc2=cos(wc2_bs);
      b_notch2=[1, -2*cos_wc2, 1];
      a_notch2=[1, -2*r*cos_wc2, r^2];

      %combinar dnv 2 filtros
      b_bs_total=conv(b_notch1, b_notch2);
      a_bs_total=conv(a_notch1, a_notch2);

      [H_bs, f_bs]=freqz(b_bs_total, a_bs_total, 1024, fs);
      H_dc=sum(b_bs_total)/sum(a_bs_total);
      b_bs_total=b_bs_total/H_dc;

      [H_bs, f_bs]=freqz(b_bs_total, a_bs_total, 1024, fs); %freq normalizada

      plot(f_bs, 20*log10(abs(H_bs)), colors(idx), 'LineWidth', 1.5,'DisplayName', ['r=', num2str(r)]);
    endfor

    hold off;
    grid on;
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title(['Filtro Rejeita-Faixa: fc1=', num2str(fc1_bs), ' Hz, fc2=', num2str(fc2_bs), ' Hz']);
    legend('show');
    xlim([0, fs/2]);
    ylim([-50, 5]);

    hold on;
    plot([fc1_bs, fc1_bs], ylim, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    plot([fc2_bs, fc2_bs], ylim, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    hold off;

  case 4
    fc_low=1000;
    fc_high=2000;
    fc1_bp=500;
    fc2_bp=2000;
    fc1_bs=1000;
    fc2_bs=2500;

    %frequencias normalizadas
    Wn_low=fc_low/(fs/2);
    Wn_high=fc_high/(fs/2);
    Wn_bp=[fc1_bp, fc2_bp]/(fs/2);
    Wn_bs=[fc1_bs, fc2_bs]/(fs/2);

    M_values=[20, 50, 100];

    windows={'rectwin', 'triang', 'bartlett', 'hamming', 'hann', 'blackman'};
    window_names={'Retangular', 'Triangular', 'Bartlett', 'Hamming', 'Hann', 'Blackman'};

    %a) passa baixa fir
    figure('Name', 'FIR Passa-Baixas - Comparacao de Janelas', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    for win_idx=1:length(windows)
      subplot(3, 2, win_idx);
      hold on;

      for M_idx=1:length(M_values)
        M=M_values(M_idx);
        win=feval(windows{win_idx}, M+1);
        b_fir=fir1(M, Wn_low, 'low', win);
        [H_fir, f_fir]=freqz(b_fir, 1, 1024, fs);

        plot(f_fir, 20*log10(abs(H_fir)), 'LineWidth', 1.5,'DisplayName', ['M=', num2str(M)]);
      endfor

      hold off;
      grid on;
      title(['Passa-Baixas - ', window_names{win_idx}]);
      xlabel('Frequencia (Hz)');
      ylabel('Magnitude (dB)');
      legend('Location', 'northeast');
      xlim([0, fs/2]);
      ylim([-100, 5]);

      hold on;
      plot([fc_low, fc_low], ylim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
      hold off;
    endfor

    %b) passa alta fir
    figure('Name', 'FIR Passa-Altas - Comparacao de Janelas', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    for win_idx=1:length(windows)
      subplot(3, 2, win_idx);
      hold on;

      for M_idx=1:length(M_values)
        M=M_values(M_idx);
        win=feval(windows{win_idx}, M+1);
        b_fir=fir1(M, Wn_high, 'high', win);
        [H_fir, f_fir]=freqz(b_fir, 1, 1024, fs);

        plot(f_fir, 20*log10(abs(H_fir)), 'LineWidth', 1.5,'DisplayName', ['M=', num2str(M)]);
      endfor

      hold off;
      grid on;
      title(['Passa-Altas - ', window_names{win_idx}]);
      xlabel('Frequencia (Hz)');
      ylabel('Magnitude (dB)');
      legend('Location', 'northeast');
      xlim([0, fs/2]);
      ylim([-100, 5]);

      hold on;
      plot([fc_high, fc_high], ylim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
      hold off;
    endfor

    %c) passa faixa fir
    figure('Name', 'FIR Passa-Faixa - Comparacao de Janelas', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    for win_idx=1:length(windows)
      subplot(3, 2, win_idx);
      hold on;

      for M_idx=1:length(M_values)
        M=M_values(M_idx);
        win=feval(windows{win_idx}, M+1);
        b_fir=fir1(M, Wn_bp, 'bandpass', win);
        [H_fir, f_fir]=freqz(b_fir, 1, 1024, fs);

        plot(f_fir, 20*log10(abs(H_fir)), 'LineWidth', 1.5,'DisplayName', ['M=', num2str(M)]);
      endfor

      hold off;
      grid on;
      title(['Passa-Faixa - ', window_names{win_idx}]);
      xlabel('Frequencia (Hz)');
      ylabel('Magnitude (dB)');
      legend('Location', 'northeast');
      xlim([0, fs/2]);
      ylim([-100, 5]);

      hold on;
      plot([fc1_bp, fc1_bp], ylim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
      plot([fc2_bp, fc2_bp], ylim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
      hold off;
    endfor

    %d) rejeita faixa fir
    figure('Name', 'FIR Rejeita-Faixa - Comparacao de Janelas', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    for win_idx=1:length(windows)
      subplot(3, 2, win_idx);
      hold on;

      for M_idx=1:length(M_values)
        M=M_values(M_idx);
        win=feval(windows{win_idx}, M+1);
        b_fir=fir1(M, Wn_bs, 'stop', win);
        [H_fir, f_fir]=freqz(b_fir, 1, 1024, fs);

        plot(f_fir, 20*log10(abs(H_fir)), 'LineWidth', 1.5,'DisplayName', ['M=', num2str(M)]);
      endfor

      hold off;
      grid on;
      title(['Rejeita-Faixa - ', window_names{win_idx}]);
      xlabel('Frequencia (Hz)');
      ylabel('Magnitude (dB)');
      legend('Location', 'northeast');
      xlim([0, fs/2]);
      ylim([-100, 5]);

      hold on;
      plot([fc1_bs, fc1_bs], ylim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
      plot([fc2_bs, fc2_bs], ylim, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
      hold off;
    endfor

  case 5
    %NOTA: acabou a luz e perdi a 5 inteira. Nice! >:]
    load('handel.mat');
    fs_handel=Fs;
    x=y;
    N=length(x);
    t=(0:N-1)/fs_handel;
    interferente1=0.05*cos(200*pi*t);
    interferente2=0.075*sin(4000*pi*t);
    sigma2=10^(-2);
    sigma=sqrt(sigma2);
    n=sigma*randn(size(x));
    y_sinal=x + interferente1' + interferente2' + n;

    Nfft=2^nextpow2(N);
    X=fft(x, Nfft);
    f_x=fs_handel*(0:Nfft/2)/Nfft;
    X_mag=abs(X(1:Nfft/2+1));
    X_mag_dB=20*log10(X_mag/max(X_mag));
    Y=fft(y_sinal, Nfft);
    Y_mag=abs(Y(1:Nfft/2+1));
    Y_mag_dB=20*log10(Y_mag/max(Y_mag));

    %NOTA: O RASCUNHO SALVA, MEU DEUS DO CÉU!

    figure('Name', 'Sinais no Dominio do Tempo', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);
    subplot(2,1,1);
    plot(t, x, 'b', 'LineWidth', 1);
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    title('Sinal Original x(t) - Handel');
    grid on;
    xlim([0, min(2, max(t))]);

    subplot(2,1,2);
    plot(t, y_sinal, 'r', 'LineWidth', 1);
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    title('Sinal Corrompido y(t) = x(t) + interferentes + ruído');
    grid on;
    xlim([0, min(2, max(t))]);

    figure('Name', 'Espectros dos Sinais', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    subplot(3,1,1);
    plot(f_x, X_mag_dB, 'b', 'LineWidth', 1.5);
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Espectro do Sinal Original x(t)');
    grid on;
    xlim([0, fs_handel/2]);
    ylim([-100, 0]);

    subplot(3,1,2);
    plot(f_x, Y_mag_dB, 'r', 'LineWidth', 1.5);
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Espectro do Sinal Corrompido y(t)');
    grid on;
    xlim([0, fs_handel/2]);
    ylim([-100, 0]);

    subplot(3,1,3);
    diff_dB=Y_mag_dB - X_mag_dB;
    plot(f_x, diff_dB, 'g', 'LineWidth', 1.5);
    xlabel('Frequencia (Hz)');
    ylabel('Diferenca (dB)');
    title('Diferenca entre Espectros: Y(f) - X(f)');
    grid on;
    xlim([0, fs_handel/2]);

    figure('Name', 'Detalhes das Interferencias', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 400]);

    f1=100;
    f2=2000;
    [~, idx1]=min(abs(f_x - f1));
    [~, idx2]=min(abs(f_x - f2));

    %NOTA: aqui to indo no "instinto", pq eu realmente n sei se isso ta certo kkkkkkkkk
    subplot(1,2,1);
    plot(f_x(max(1,idx1-50):min(length(f_x),idx1+50)),Y_mag_dB(max(1,idx1-50):min(length(f_x),idx1+50)), 'r', 'LineWidth', 2);
    hold on;
    plot(f_x(max(1,idx1-50):min(length(f_x),idx1+50)),X_mag_dB(max(1,idx1-50):min(length(f_x),idx1+50)), 'b--', 'LineWidth', 1.5);
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title(['Detalhe em f = ', num2str(f1), ' Hz (interferência cos)']);
    grid on;
    legend('y(t)', 'x(t)');

    subplot(1,2,2);
    plot(f_x(max(1,idx2-50):min(length(f_x),idx2+50)),Y_mag_dB(max(1,idx2-50):min(length(f_x),idx2+50)), 'r', 'LineWidth', 2);
    hold on;
    plot(f_x(max(1,idx2-50):min(length(f_x),idx2+50)),X_mag_dB(max(1,idx2-50):min(length(f_x),idx2+50)), 'b--', 'LineWidth', 1.5);
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title(['Detalhe em f = ', num2str(f2), ' Hz (interferencia sin)']);
    grid on;
    legend('y(t)', 'x(t)');

    % ------------------------- QUESTÃO 6 ---------------------------------

    f_notch1=100;
    f_notch2=2000;
    fc_lp=3500;

    Wn_notch1=f_notch1/(fs_handel/2);
    Wn_notch2=f_notch2/(fs_handel/2);
    Wn_lp=fc_lp/(fs_handel/2);

    r_notch=0.98;
    wc_notch1=2*pi*f_notch1/fs_handel;
    cos_wc1=cos(wc_notch1);
    b_notch1=[1, -2*cos_wc1, 1];
    a_notch1=[1, -2*r_notch*cos_wc1, r_notch^2];
    wc_notch2=2*pi*f_notch2/fs_handel;
    cos_wc2=cos(wc_notch2);
    b_notch2=[1, -2*cos_wc2, 1];
    a_notch2=[1, -2*r_notch*cos_wc2, r_notch^2];

    M_fir=100;
    win=hamming(M_fir+1);
    b_lp=fir1(M_fir, Wn_lp, 'low', win);
    a_lp=1;

    %combina notches
    b_notches=conv(b_notch1, b_notch2);
    a_notches=conv(a_notch1, a_notch2);

    %combinar com passa baixa
    b_total=conv(b_notches, b_lp);
    a_total=conv(a_notches, a_lp);

    [H_test, ~]=freqz(b_total, a_total, 1024, fs_handel);
    H_max=max(abs(H_test));
    b_total=b_total/H_max;

    figure('Name', 'Questao 6 - Resposta do Filtro Projetado', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);

    [H_filter, f_filter]=freqz(b_total, a_total, 1024, fs_handel);

    subplot(2,2,1);
    plot(f_filter, 20*log10(abs(H_filter)), 'b', 'LineWidth', 2);
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Resposta em Magnitude do Filtro');
    grid on;
    xlim([0, fs_handel/2]);
    ylim([-80, 5]);

    hold on;
    plot([f_notch1, f_notch1], ylim, 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot([f_notch2, f_notch2], ylim, 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');
    plot([fc_lp, fc_lp], ylim, 'g--', 'LineWidth', 1, 'HandleVisibility', 'off');
    hold off;

    subplot(2,2,2);
    plot(f_filter, unwrap(angle(H_filter)), 'b', 'LineWidth', 2);
    xlabel('Frequencia (Hz)');
    ylabel('Fase (radianos)');
    title('Resposta em Fase do Filtro');
    grid on;
    xlim([0, fs_handel/2]);

    subplot(2,2,3);
    zplane(b_total, a_total);
    title('Diagrama de Polos e Zeros');
    grid on;

    subplot(2,2,4);
    imp_resp=filter(b_total, a_total, [1; zeros(100,1)]);
    stem(0:100, imp_resp, 'b', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlabel('Amostras (n)');
    ylabel('Amplitude');
    title('Resposta ao Impulso (100 primeiras amostras)');
    grid on;

    x_recuperado=filter(b_total, a_total, y_sinal);

    figure('Name', 'Questao 6 - Comparacaoo de Espectros', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 800]);

    X_rec=fft(x_recuperado, Nfft);
    X_rec_mag=abs(X_rec(1:Nfft/2+1));
    X_rec_mag_dB=20*log10(X_rec_mag/max(X_rec_mag));

    subplot(3,1,1);
    plot(f_x, X_mag_dB, 'b', 'LineWidth', 1.5);
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Espectro do Sinal Original x(t)');
    grid on;
    xlim([0, fs_handel/2]);
    ylim([-100, 0]);

    subplot(3,1,2);
    plot(f_x, Y_mag_dB, 'r', 'LineWidth', 1.5);
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Espectro do Sinal Corrompido y(t)');
    grid on;
    xlim([0, fs_handel/2]);
    ylim([-100, 0]);

    subplot(3,1,3);
    plot(f_x, X_rec_mag_dB, 'g', 'LineWidth', 1.5);
    xlabel('Frequencia (Hz)');
    ylabel('Magnitude (dB)');
    title('Espectro do Sinal Recuperado apos Filtragem');
    grid on;
    xlim([0, fs_handel/2]);
    ylim([-100, 0]);

    figure('Name', 'Questao 6 - Sinais no Tempo', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 600]);

    t_seg=0.5;
    n_amostras=min(round(t_seg*fs_handel), N);

    subplot(3,1,1);
    plot(t(1:n_amostras), x(1:n_amostras), 'b', 'LineWidth', 1.5);
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    title('Sinal Original x(t)');
    grid on;
    xlim([0, t_seg]);

    subplot(3,1,2);
    plot(t(1:n_amostras), y_sinal(1:n_amostras), 'r', 'LineWidth', 1.5);
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    title('Sinal Corrompido y(t)');
    grid on;
    xlim([0, t_seg]);

    subplot(3,1,3);
    plot(t(1:n_amostras), x_recuperado(1:n_amostras), 'g', 'LineWidth', 1.5);
    xlabel('Tempo (s)');
    ylabel('Amplitude');
    title('Sinal Recuperado após Filtragem');
    grid on;
    xlim([0, t_seg]);

    try
      disp('Sinal original x(t)');
      pause(1);
      soundsc(x, fs_handel);
      pause(2);

      disp('Sinal corrompido y(t)'); %NOTA: lampada de fluorescente
      pause(1);
      soundsc(y_sinal, fs_handel);
      pause(2);

      disp('Sinal recuperado apos filtragem'); %NOta: internet discada
      pause(1);
      soundsc(x_recuperado, fs_handel);
      pause(2);
    catch
      disp('Reproducao de audio nao disponivel.');
    end_try_catch

  otherwise
    disp("Opcao invalida. Terminando o Programa.");
    pause(2);
end
