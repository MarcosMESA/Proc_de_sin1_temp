%Aluno: Marcos Mauricio da Cunha de Morais

clear
clc
close all;


try
    pkg load signal;
    load handel.mat
    %load("handel", "-mat");
catch
    disp("Arquivo não encontrado.");
    pkg install -forge signal
    pkg load signal
end

try
  disp("Carregando audio. ");
  pkg load audio
catch
  disp("Audio nao encontrado. Instalando audio. ");
  pkg install -forge audio
  pkg load audio
end

prompt='Qual Exercicio quer simular? Digite o Numero de 1 a 5: ';
option=input(prompt);

switch option
    case 1
      taxas=[0.995, 0.99, 0.90, 0.75, 0.50]
      num_taxas=length(taxas)
      Y_dct=dct(y)
      N=length(Y_dct)

      eqm_values=zeros(1, num_taxas);
      num_coefs=zeros(1, num_taxas);
      taxas_compressao=zeros(1, num_taxas);

      figure('Position',[100, 100, 1200, 800]);

      for k=1:num_taxas
        taxa=taxas(k);

        %ordenar coeficiente por magnitude
        [coef_ordenados,indices]=sort(abs(Y_dct),'descend');

        %num de coef necessarios
        energia_total=norm(Y_dct)^2;
        energia_acumulada=0;
        num_coef=0;

        while (energia_acumulada / energia_total) < taxa && num_coef < N
          num_coef=num_coef+1;
          energia_acumulada=energia_acumulada+abs(Y_dct(indices(num_coef)))^2;
        end

        Y_comprimido=Y_dct;
        Y_comprimido(indices(num_coef+1:end))=0;

        %reconstroi sinal
        y_comprimido=idct(Y_comprimido);

        %calcula estatisticas
        taxa_compressao_val=(1 - (num_coef / N))*100;

        %eqm
        eqm=mean((y - y_comprimido).^2);
        erro_quadratico=norm(y - y_comprimido)^2 / norm(y)^2;

        eqm_values(k) = eqm;
        erro_quadratico_values(k) = erro_quadratico;
        num_coefs(k) = num_coef;
        taxas_compressao(k) = taxa_compressao_val;

        %plota
        subplot(num_taxas,2,(k-1)*2 + 1);
        plot(y,'b','LineWidth',1);
        hold on;
        plot(y_comprimido,'r--','LineWidth',1);
        title(sprintf('Sinal - %.1f%% Compressao',taxa*100));
        xlabel('Amostras');
        ylabel('Amplitude');
        legend('Original','Comprimido','Location','northeast');
        grid on;

        subplot(num_taxas,2,(k-1)*2 + 2);
        stem(1:N,abs(Y_dct),'b.', 'MarkerSize',1);
        hold on;
        stem(indices(1:num_coef),abs(Y_dct(indices(1:num_coef))),'r.', 'MarkerSize',2);
        title(sprintf('Coeficientes DCT - %d usados',num_coef));
        xlabel('Indice do Coeficiente');
        ylabel('Magnitude');
        legend('Todos Coef.','Coef. Mantidos','Location','northeast');
        grid on;

        fprintf('Taxa de Compressão: %.1f%% ', taxa*100);
        fprintf('Coeficientes usados: %d de %d\n', num_coef, N);
        fprintf('Taxa de compressão: %.2f%%\n', taxa_compressao_val);
        fprintf('Erro Quadrático Médio (EQM): %.6f\n', eqm);
        fprintf('Erro Quadrático Relativo: %.6f\n', erro_quadratico);
      end

      sound(y, Fs);
      pause(length(y)/Fs + 1);

      for k=1:num_taxas
        %recalcular sinal comprimido pra reproduzie
        Y_comprimido=Y_dct;
        [~,indices]=sort(abs(Y_dct),'descend');

        energia_total=norm(Y_dct)^2;
        energia_acumulada=0;
        num_coef=0;

        while (energia_acumulada/energia_total)<taxas(k) && num_coef<N
          num_coef=num_coef + 1;
          energia_acumulada=energia_acumulada+abs(Y_dct(indices(num_coef)))^2;
        end

        Y_comprimido(indices(num_coef+1:end))=0;
        y_comprimido=idct(Y_comprimido);

        sound(y_comprimido, Fs);
        pause(length(y_comprimido)/Fs+1);
      end

      %NOTA: esse aqui foi 3 tiros na cabeça

    case 2
      taxas=[0.995, 0.99, 0.90, 0.75, 0.50];
      num_taxas=length(taxas);
      Y_fft=fft(y);
      N=length(Y_fft);

      eqm_values=zeros(1, num_taxas);
      erro_quadratico_values=zeros(1, num_taxas);
      num_coefs_total=zeros(1, num_taxas);
      taxas_compressao=zeros(1, num_taxas);

      figure('Position', [100, 100, 1200, 800]);

      for k=1:num_taxas
        taxa=taxas(k);
        magnitudes=abs(Y_fft(1:floor(N/2)+1));
        [coef_ordenados,indices]=sort(magnitudes,'descend');
        energia_total=norm(Y_fft)^2;
        energia_acumulada=0;
        num_coef=0;

        %considera coeficiente simetricos pra simetria da fft
        while (energia_acumulada/energia_total)<taxa && num_coef<floor(N/2)
          num_coef=num_coef+1;
          idx=indices(num_coef);
          %adiciona energia do coeficiente e seu simetrico
          if idx==1
            energia_acumulada=energia_acumulada+2*abs(Y_fft(idx))^2;
          else
            energia_acumulada=energia_acumulada+2*abs(Y_fft(idx))^2;
          end
        end
        Y_comprimido=zeros(size(Y_fft));

        for j=1:num_coef
          idx=indices(j);
          Y_comprimido(idx)=Y_fft(idx);
          if idx>1 && idx<=floor(N/2)
            %coeficiente simetrico
            idx_simetrico=N-idx+2;
            Y_comprimido(idx_simetrico)=Y_fft(idx_simetrico);
          end
        end

        y_comprimido=real(ifft(Y_comprimido));

        num_coef_total=2*num_coef - 1;
        if num_coef_total>N
          num_coef_total=N;
        end

        taxa_compressao_val=(1-(num_coef_total/N))*100;
        taxa_compressao=(1-(num_coef_total/N))*100;

        eqm = mean((y - y_comprimido).^2);
        erro_quadratico=norm(y - y_comprimido)^2/norm(y)^2;

        eqm_values(k)=eqm;
        erro_quadratico_values(k)=erro_quadratico;
        num_coefs_total(k)=num_coef_total;
        taxas_compressao(k)=taxa_compressao_val;

        subplot(num_taxas,2,(k-1)*2+1);
        plot(y,'b','LineWidth',1);
        hold on;
        plot(y_comprimido,'r--','LineWidth', 1);
        title(sprintf('Sinal DFT - %.1f%% Compressao\nEQM: %.4f',taxa*100,eqm));
        xlabel('Amostras');
        ylabel('Amplitude');
        legend('Original','Comprimido','Location', 'northeast');
        grid on;

        subplot(num_taxas,2,(k-1)*2 + 2);
        freq=(0:N-1)*Fs/N;
        stem(freq(1:floor(N/2)),abs(Y_fft(1:floor(N/2))),'b.','MarkerSize',1);
        hold on;

        coef_mantidos = zeros(1, floor(N/2));
        for j=1:num_coef
          idx=indices(j);
          if idx<=floor(N/2)
            coef_mantidos(idx)=abs(Y_fft(idx));
          end
        end

        stem(freq(1:floor(N/2)),coef_mantidos, 'r.','MarkerSize',2);

        title(sprintf('Coeficientes FFT - %d pares usados',num_coef));
        xlabel('Frequencia (Hz)');
        ylabel('Magnitude');
        legend('Todos Coef.','Coef. Mantidos','Location','northeast');
        grid on;

        fprintf('\n--- DFT - Taxa de Compressão: %.1f%% ---\n', taxa*100);
        fprintf('Coeficientes usados: %d de %d\n', num_coef_total, N);
        fprintf('Taxa de compressão: %.2f%%\n', taxa_compressao_val);
        fprintf('Erro Quadratico Médio (EQM): %.6f\n', eqm);
        fprintf('Erro Quadratico Relativo: %.6f\n', erro_quadratico);
      end

      figure('Position', [150, 150, 800, 400]);
      subplot(1, 2, 1);
      plot(taxas*100, eqm_values, 'ro-', 'LineWidth', 2, 'MarkerSize', 8);
      xlabel('Taxa de Compressao Alvo (%)');
      ylabel('Erro Quadratico Medio (EQM)');
      title('DFT - EQM vs Taxa de Compressao');
      grid on;

      subplot(1, 2, 2);
      semilogy(taxas*100, eqm_values, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
      xlabel('Taxa de Compressao Alvo (%)');
      ylabel('EQM (escala log)');
      title('DFT - EQM vs Taxa de Compressao (Escala Log)');
      grid on;

      figura_comparacao=figure('Position',[100, 100, 1400, 600]);

      taxa_comparacao=0.50;

      Y_dct=dct(y);
      [~, indices_dct]=sort(abs(Y_dct),'descend');

      energia_total_dct=norm(Y_dct)^2;
      energia_acum_dct=0;
      num_coef_dct=0;

      while (energia_acum_dct/energia_total_dct)<taxa_comparacao && num_coef_dct<N
        num_coef_dct=num_coef_dct+1;
        energia_acum_dct=energia_acum_dct+abs(Y_dct(indices_dct(num_coef_dct)))^2;
      end

      Y_dct_comprimido=Y_dct;
      Y_dct_comprimido(indices_dct(num_coef_dct+1:end))=0;
      y_dct_comprimido=idct(Y_dct_comprimido);

      [~, indices_fft]=sort(abs(Y_fft(1:floor(N/2)+1)),'descend');
      energia_total_fft=norm(Y_fft)^2;
      energia_acum_fft=0;
      num_coef_fft=0;

      while (energia_acum_fft/energia_total_fft)<taxa_comparacao && num_coef_fft<floor(N/2)
        num_coef_fft=num_coef_fft+1;
        idx = indices_fft(num_coef_fft);
        if idx==1
          energia_acum_fft=energia_acum_fft+2*abs(Y_fft(idx))^2;
        else
          energia_acum_fft=energia_acum_fft+2*abs(Y_fft(idx))^2;
        end
      end

      Y_fft_comprimido=zeros(size(Y_fft));
      for j=1:num_coef_fft
        idx=indices_fft(j);
        Y_fft_comprimido(idx)=Y_fft(idx);
        if idx>1 && idx<=floor(N/2)
          idx_simetrico=N - idx+2;
          Y_fft_comprimido(idx_simetrico)=Y_fft(idx_simetrico);
        end
      end
      y_fft_comprimido=real(ifft(Y_fft_comprimido));

      % Cálculo dos EQM para comparação - ADICIONADO
      eqm_dct = mean((y - y_dct_comprimido).^2);
      eqm_fft = mean((y - y_fft_comprimido).^2);

      %plot
      subplot(2, 3, 1);
      plot(y, 'k', 'LineWidth', 2);
      title('Sinal Original');
      xlabel('Amostras');
      ylabel('Amplitude');
      grid on;

      subplot(2, 3, 2);
      plot(y_dct_comprimido, 'b', 'LineWidth', 1.5);
      title(sprintf('DCT Comprimido (%d coef)\nEQM: %.4f',num_coef_dct,eqm_dct));
      xlabel('Amostras');
      ylabel('Amplitude');
      grid on;

      subplot(2, 3, 3);
      plot(y_fft_comprimido, 'r', 'LineWidth', 1.5);
      title(sprintf('DFT Comprimido (%d coef)\nEQM: %.4f',2*num_coef_fft-1,eqm_fft));
      xlabel('Amostras');
      ylabel('Amplitude');
      grid on;

      subplot(2, 3, 4);
      plot(y - y_dct_comprimido, 'b', 'LineWidth', 1);
      title('Erro DCT');
      xlabel('Amostras');
      ylabel('Amplitude');
      grid on;

      subplot(2, 3, 5);
      plot(y - y_fft_comprimido, 'r', 'LineWidth', 1);
      title('Erro DFT');
      xlabel('Amostras');
      ylabel('Amplitude');
      grid on;

      subplot(2, 3, 6);
      erro_dct=norm(y - y_dct_comprimido)^2 /norm(y)^2;
      erro_fft=norm(y - y_fft_comprimido)^2 /norm(y)^2;
       bar([eqm_dct, eqm_fft]);
      set(gca, 'XTickLabel', {'DCT', 'DFT'});
      title('Erro Quadratico Medio (EQM)');
      ylabel('EQM');
      grid on;

      fprintf('COMPARAÇAO DCT vs DFT (50%% compressao)');
      fprintf('DCT - Coeficientes: %d, EQM: %.6f\n',num_coef_dct,eqm_dct);
      fprintf('DFT - Coeficientes: %d, EQM: %.6f\n',2*num_coef_fft-1,eqm_fft);

      sound(y,Fs);
      pause(length(y)/Fs+1);
      sound(y_dct_comprimido,Fs);
      pause(length(y_dct_comprimido)/Fs+1);
      sound(y_fft_comprimido,Fs);
      pause(length(y_fft_comprimido)/Fs+1);

      %NOTA: Es hora de mimir

    case 3
      try
        imagem=imread('sosias.jpg');
        %imagem=imread('robot_masters_rockman1_2.jpg');

        [altura, largura, canais] = size(imagem);

        if canais==3
          imagem_cinza=rgb2gray(imagem);
          imagem_original=imagem_cinza;
        else
          imagem_original=imagem;
        end

        imagem_double=double(imagem_original);
        imagem_dct=dct2(imagem_double);

        %calcular energia da imagem
        energia_total = sum(imagem_double(:).^2);
        energia_dct = sum(imagem_dct(:).^2);

        fprintf('Energia da imagem original: %.4f\n', energia_total);
        fprintf('Energia da imagem DCT: %.4f\n', energia_dct);
        fprintf('Diferença de energia: %.6f\n', abs(energia_total - energia_dct));

        figure('Position', [100, 100, 1000, 400]);

        subplot(1, 3, 1);
        imshow(imagem);
        title('Imagem Original');

        subplot(1, 3, 2);
        imshow(imagem_original);
        title('Imagem Original em Preto e Branco');

        subplot(1, 3, 3);
        imshow(log(abs(imagem_dct) + 1), []);
        title('DCT 2D (escala log)');
        colorbar;

      catch
        fprintf('Imagem nao encontrada');
      end

    case 4
      try
        imagem=imread('sosias.jpg');
        %imagem=imread('robot_masters_rockman1_2.jpg');

        [altura, largura, canais]=size(imagem);

        if canais==3
          imagem_cinza=rgb2gray(imagem);
          imagem_original=imagem_cinza;
        else
          imagem_original=imagem;
        end

        imagem_double=double(imagem_original);
        imagem_dct=dct2(imagem_double);

        energia_total=sum(imagem_double(:).^2);
        energia_dct_total=sum(imagem_dct(:).^2);

        fprintf('Dimensoes da imagem: %d x %d pixels\n', altura, largura);
        fprintf('Total de coeficientes DCT: %d\n', numel(imagem_dct));
        fprintf('Energia total da imagem: %.4f\n', energia_total);

        taxas=[0.995, 0.99, 0.90, 0.75, 0.50];
        num_taxas=length(taxas);

        num_coefs_mantidos=zeros(1, num_taxas);
        taxas_compressao=zeros(1, num_taxas);
        eqm_values=zeros(1, num_taxas);

        figure('Position',[100, 100, 1400, 800]);

        for k=1:num_taxas
          taxa_energia=taxas(k);

          coef_ordenados=sort(abs(imagem_dct(:)), 'descend');
          energia_acumulada=cumsum(coef_ordenados.^2);

          num_coef=find(energia_acumulada >= taxa_energia * energia_dct_total, 1);

          if isempty(num_coef)
            num_coef=numel(imagem_dct);
          end

          limiar=coef_ordenados(num_coef);
          mascara=abs(imagem_dct) >= limiar;

          imagem_dct_comprimida=imagem_dct .* mascara;

          imagem_reconstruida=idct2(imagem_dct_comprimida);
          imagem_reconstruida=uint8(max(0, min(255, imagem_reconstruida)));

          num_coefs_mantidos(k)=sum(mascara(:));
          taxas_compressao(k)=(1 - num_coefs_mantidos(k) / numel(imagem_dct)) * 100;
          eqm_values(k)=mean((double(imagem_original(:)) - double(imagem_reconstruida(:))).^2);

          subplot(3, num_taxas, k);
          imshow(imagem_reconstruida);
          title(sprintf('%.1f%% Energia\n%d coef. (%.1f%%)', ...
          taxa_energia*100, num_coefs_mantidos(k), taxas_compressao(k)));

          subplot(3, num_taxas, num_taxas + k);
          mascara_visual=double(mascara);
          imshow(mascara_visual);
          title(sprintf('Coef. Mantidos: %d/%d', num_coefs_mantidos(k), numel(imagem_dct)));

          subplot(3, num_taxas, 2*num_taxas + k);
          diferenca=double(imagem_original) - double(imagem_reconstruida);
          imshow(diferenca, []);
          title(sprintf('Erro (EQM: %.2f)', eqm_values(k)));
          colorbar;

          fprintf('Compressao para %.1f%% da energia', taxa_energia*100);
          fprintf('Coeficientes mantidos: %d de %d\n', num_coefs_mantidos(k), numel(imagem_dct));
          fprintf('Taxa de compressao: %.2f%%\n', taxas_compressao(k));
          fprintf('Erro Quadratico Medio: %.4f\n', eqm_values(k));
        end

        figure('Position',[150, 150, 1000, 600]);

        subplot(2,2,1);
        coef_ordenados=sort(abs(imagem_dct(:)), 'descend');
        energia_acumulada=cumsum(coef_ordenados.^2) / energia_dct_total * 100;
        plot(energia_acumulada, 'b-', 'LineWidth',2);
        xlabel('Numero de Coeficientes');
        ylabel('Energia Acumulada (%)');
        title('Energia vs Coeficientes DCT');
        grid on;

        hold on;
        for k=1:num_taxas
          plot([num_coefs_mantidos(k), num_coefs_mantidos(k)], [0, taxas(k)*100], 'r--');
          plot([1, num_coefs_mantidos(k)], [taxas(k)*100, taxas(k)*100], 'r--');
        end

        subplot(2,2,2);
        plot(taxas*100, num_coefs_mantidos, 'ro-', 'LineWidth',2, 'MarkerSize',8);
        xlabel('Energia Mantida (%)');
        ylabel('Coeficientes Necessarios');
        title('Coeficientes vs Energia');
        grid on;

        subplot(2,2,3);
        plot(taxas*100, taxas_compressao, 'go-', 'LineWidth',2, 'MarkerSize',8);
        xlabel('Energia Mantida (%)');
        ylabel('Taxa de Compressao (%)');
        title('Compressao vs Energia');
        grid on;

        subplot(2,2,4);
        semilogy(taxas*100, eqm_values, 'bo-', 'LineWidth',2, 'MarkerSize',8);
        xlabel('Energia Mantida (%)');
        ylabel('EQM (escala log)');
        title('Erro vs Energia');
        grid on;

      catch
        fprintf('Imagem nao encontrada');
      end

    case 5
      try
        imagem=imread('sosias.jpg');
        %imagem=imread('robot_masters_rockman1_2.jpg');

        [altura, largura, canais]=size(imagem);

        if canais==3
          imagem_cinza=rgb2gray(imagem);
          imagem_original=imagem_cinza;
        else
          imagem_original=imagem;
        end

        imagem_double=double(imagem_original);

        blocos_tamanhos=[8, 16];
        num_tamanhos=length(blocos_tamanhos);

        taxas=[0.995, 0.99, 0.90, 0.75, 0.50];
        num_taxas=length(taxas);

        figure('Position',[100, 100, 1400, 1000]);

        for bloco_idx=1:num_tamanhos
          tamanho_bloco=blocos_tamanhos(bloco_idx);
          fprintf('Processando blocos de %dx%d', tamanho_bloco, tamanho_bloco);

          num_blocos_vertical=floor(altura/tamanho_bloco);
          num_blocos_horizontal=floor(largura/tamanho_bloco);

          fprintf('Numero de blocos: %d x %d = %d blocos\n',num_blocos_vertical, num_blocos_horizontal,num_blocos_vertical*num_blocos_horizontal);

          imagem_reconstruida=zeros(size(imagem_double));

          for i=1:num_blocos_vertical
            for j=1:num_blocos_horizontal
              linha_inicio=(i-1)*tamanho_bloco+1;
              linha_fim=i*tamanho_bloco;
              coluna_inicio=(j-1)*tamanho_bloco+1;
              coluna_fim=j*tamanho_bloco;

              bloco=imagem_double(linha_inicio:linha_fim, coluna_inicio:coluna_fim);

              bloco_dct=dct2(bloco);

              energia_bloco_total=sum(bloco(:).^2);
              energia_bloco_dct=sum(bloco_dct(:).^2);

              coef_ordenados=sort(abs(bloco_dct(:)), 'descend');
              energia_acumulada=cumsum(coef_ordenados.^2);

              taxa_energia=0.90;
              num_coef=find(energia_acumulada >= taxa_energia * energia_bloco_dct, 1);

              if isempty(num_coef)
                num_coef=numel(bloco_dct);
              end

              limiar=coef_ordenados(num_coef);
              mascara=abs(bloco_dct) >= limiar;
              bloco_dct_comprimido=bloco_dct .* mascara;

              bloco_reconstruido=idct2(bloco_dct_comprimido);

              imagem_reconstruida(linha_inicio:linha_fim, coluna_inicio:coluna_fim)=bloco_reconstruido;
            end
          end

          imagem_reconstruida_uint8=uint8(max(0, min(255, imagem_reconstruida)));
          eqm=mean((imagem_double(:) - imagem_reconstruida(:)).^2);
          psnr=10*log10(255^2/eqm);

          subplot(3, num_tamanhos, bloco_idx);
          imshow(imagem_reconstruida_uint8);
          title(sprintf('Blocos %dx%d\nEQM: %.2f, PSNR: %.2f dB',tamanho_bloco, tamanho_bloco, eqm, psnr));

          subplot(3, num_tamanhos, num_tamanhos + bloco_idx);
          bloco_exemplo=imagem_double(1:tamanho_bloco, 1:tamanho_bloco);
          bloco_exemplo_dct=dct2(bloco_exemplo);
          imshow(log(abs(bloco_exemplo_dct)+1), []);
          title(sprintf('DCT Bloco %dx%d Exemplo', tamanho_bloco, tamanho_bloco));
          colorbar;

          subplot(3, num_tamanhos, 2*num_tamanhos + bloco_idx);
          diferenca=imagem_double - imagem_reconstruida;
          imshow(diferenca, []);
          title(sprintf('Erro - Blocos %dx%d', tamanho_bloco, tamanho_bloco));
          colorbar;

          fprintf('EQM: %.4f\n', eqm);
          fprintf('PSNR: %.2f dB\n', psnr);
        end

        imagem_dct_global=dct2(imagem_double);
        coef_ordenados_global=sort(abs(imagem_dct_global(:)), 'descend');
        energia_total_global=sum(imagem_dct_global(:).^2);
        energia_acumulada_global=cumsum(coef_ordenados_global.^2);

        taxa_comparacao=0.90;
        num_coef_global=find(energia_acumulada_global >= taxa_comparacao * energia_total_global, 1);
        limiar_global=coef_ordenados_global(num_coef_global);
        mascara_global=abs(imagem_dct_global) >= limiar_global;
        imagem_dct_global_comprimida=imagem_dct_global .* mascara_global;
        imagem_global_reconstruida=idct2(imagem_dct_global_comprimida);
        imagem_global_reconstruida_uint8=uint8(max(0, min(255, imagem_global_reconstruida)));
        eqm_global=mean((imagem_double(:) - imagem_global_reconstruida(:)).^2);
        psnr_global=10*log10(255^2/eqm_global);

        fprintf('Compressao Glibal (90%% energia):\n');
        fprintf('Coeficientes mantidos: %d de %d (%.1f%%)\n',sum(mascara_global(:)), numel(imagem_dct_global),(1-sum(mascara_global(:))/numel(imagem_dct_global))*100);
        fprintf('EQM: %.4f, PSNR: %.2f dB\n', eqm_global, psnr_global);

        figure('Position',[200, 200, 1200, 400]);

        subplot(1,3,1);
        imshow(imagem_original);
        title('Original');

        subplot(1,3,2);
        imshow(imagem_global_reconstruida_uint8);
        title(sprintf('Global\nEQM: %.2f, PSNR: %.2f dB', eqm_global, psnr_global));

        subplot(1,3,3);
        imagem_reconstruida_blocos=zeros(size(imagem_double));
        tamanho_bloco=8;
        num_blocos_vertical=floor(altura/tamanho_bloco);
        num_blocos_horizontal=floor(largura/tamanho_bloco);

        for i=1:num_blocos_vertical
          for j=1:num_blocos_horizontal
            linha_inicio=(i-1)*tamanho_bloco+1;
            linha_fim=i*tamanho_bloco;
            coluna_inicio=(j-1)*tamanho_bloco+1;
            coluna_fim=j*tamanho_bloco;

            bloco=imagem_double(linha_inicio:linha_fim, coluna_inicio:coluna_fim);
            bloco_dct=dct2(bloco);

            coef_ordenados=sort(abs(bloco_dct(:)), 'descend');
            energia_bloco_dct=sum(bloco_dct(:).^2);
            energia_acumulada=cumsum(coef_ordenados.^2);

            num_coef=find(energia_acumulada >= 0.90 * energia_bloco_dct, 1);
            limiar=coef_ordenados(num_coef);
            mascara=abs(bloco_dct) >= limiar;
            bloco_dct_comprimido=bloco_dct .* mascara;
            bloco_reconstruido=idct2(bloco_dct_comprimido);

            imagem_reconstruida_blocos(linha_inicio:linha_fim, coluna_inicio:coluna_fim)=bloco_reconstruido;
          end
        end

        imagem_reconstruida_blocos_uint8=uint8(max(0, min(255, imagem_reconstruida_blocos)));
        eqm_blocos=mean((imagem_double(:) - imagem_reconstruida_blocos(:)).^2);
        psnr_blocos=10*log10(255^2/eqm_blocos);

        imshow(imagem_reconstruida_blocos_uint8);
        title(sprintf('Blocos 8x8\nEQM: %.2f, PSNR: %.2f dB', eqm_blocos, psnr_blocos));

      catch
          fprintf('Imagem nao encontrada');
      end

      %NOTA: ACABOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOU!
      %NOTA: Pa papá papá papá
      %NOTA: Pa papá papá papá papapa
      %NOTA: Pa papá paaaaa
      %NOTA: Pa papá papá papá

    otherwise
        disp("questao nao existe. Saindo...");
end
