% Criar arquivo handel.mat para Octave
function create_handel_mat()
    % Parâmetros do sinal (baseados no handel.mat original)
    Fs = 8192;                    % Taxa de amostragem
    duration = 8.925;             % Duração em segundos (73113 amostras / 8192 Hz)
    t = 0:1/Fs:duration;
    t = t(1:end-1);               % Ajustar para ter 73113 amostras

    % Criar um sinal complexo similar ao Handel
    % Baseado em informações do sinal original
    y = zeros(1, length(t));

    % Adicionar componentes de frequência típicas de música coral
    frequencies = [196, 262, 330, 392, 494, 587, 659, 784]; % G3, C4, E4, G4, B4, D5, E5, G5
    amplitudes = [0.3, 0.8, 0.6, 1.0, 0.7, 0.5, 0.4, 0.3];

    for i = 1:length(frequencies)
        % Adicionar variação de amplitude para simular "fraseado" musical
        amp_variation = 0.7 + 0.3 * sin(2*pi*0.2*t + i*pi/4);
        y = y + amplitudes(i) * amp_variation .* sin(2*pi*frequencies(i)*t);
    end

    % Adicionar harmônicos
    y = y + 0.2 * sin(2*pi*2*392*t) + 0.1 * sin(2*pi*3*392*t);

    % Aplicar envelope de amplitude para simular início e fim de nota
    attack = 0.1; % 100ms de ataque
    release = 0.5; % 500ms de release
    envelope = ones(size(t));
    envelope(t < attack) = t(t < attack)/attack;
    envelope(t > (duration - release)) = 1 - (t(t > (duration - release)) - (duration - release))/release;

    y = y .* envelope;

    % Adicionar reverberação suave
    for i = 1:3
        delay = round(0.05*Fs*i); % Atrasos de 50ms, 100ms, 150ms
        if delay < length(y)
            y(delay+1:end) = y(delay+1:end) + 0.2^i * y(1:end-delay);
        end
    end

    % Normalizar
    y = y / max(abs(y));

    % Garantir que temos exatamente 73113 amostras (como o original)
    if length(y) > 73113
        y = y(1:73113);
    else
        y = [y, zeros(1, 73113 - length(y))];
    end

    % Salvar as variáveis no arquivo handel.mat
    save -mat handel.mat y Fs

    fprintf('Arquivo handel.mat criado com sucesso!\n');
    fprintf('Variáveis salvas: y (%d amostras) e Fs (%d Hz)\n', length(y), Fs);
end

% Executar a função para criar o arquivo
create_handel_mat();
