function synchronization_FHN_error()
    clc;
    clear;

    % Параметры моделирования
    tspan = [0, 40000];      % Время моделирования
    dt = 0.01;               % Шаг интегрирования

    % Параметры первого генератора
    etta1 = 0.5;             % Параметр etta для первого генератора
    epsilon1 = 0.0049057105; % Параметр epsilon для первого генератора

    % Параметры второго генератора
    etta2 = 0.47;            % Параметр etta для второго генератора
    epsilon2 = 0.005;        % Параметр epsilon для второго генератора

    % Начальные условия
    y10 = [0.2; 0; 0.02];    % Начальные значения для первого генератора
    y20 = [0.4; 0; 0.05];    % Начальные значения для второго генератора

    % Параметры сопротивления мемристора
    R_memristor_min_region1 = 7e2; % Минимальное сопротивление для региона 1
    R_memristor_max_region1 = 1e4; % Максимальное сопротивление для региона 1
    R_memristor_min_region2 = 2e2; % Минимальное сопротивление для региона 2
    R_memristor_max_region2 = 1.2e4; % Максимальное сопротивление для региона 2

    % Численное решение с помощью ODE45
    [t, y] = ode45(@(t, y) coupled_FHN_multi_regions(t, y, etta1, epsilon1, ...
        etta2, epsilon2, R_memristor_min_region1, R_memristor_max_region1, ...
        R_memristor_min_region2, R_memristor_max_region2), tspan, [y10; y20]);

    % Расчёт ошибки синхронизации
    error_sync = calculate_period_error(t, y);

   % Визуализация результатов в одном графике
figure;

% 1. Сигналы генераторов.
subplot(2, 1, 1);  % Первый подграфик из двух
plot(t, y(:, 1), 'b-', 'DisplayName', 'Генератор 1'); hold on;
plot(t, y(:, 4), 'r-', 'DisplayName', 'Генератор 2');
title('Сигналы генераторов');
xlabel('Время');
ylabel('Напряжение');
legend('Генератор 1', 'Генератор 2');
grid on;

% 2. Ошибка синхронизации (разность периодов).
subplot(2, 1, 2);  % Второй подграфик из двух
plot(t(1:length(error_sync)), error_sync, 'k-', 'DisplayName', 'Ошибка синхронизации');
title('Ошибка синхронизации (разность периодов)');
xlabel('Время');
ylabel('Разность периодов');
legend('Ошибка синхронизации');
grid on;

end

function dy = coupled_FHN_multi_regions(t, y, etta1, epsilon1, etta2, ...
    epsilon2, R_min1, R_max1, R_min2, R_max2)
    % Взаимодействие двух генераторов с мемристивным взаимодействием

    % Различные сопротивления мемристора в зависимости от времени
    if t < 50  % Первая область взаимодействия
        R_min = R_min1;
        R_max = R_max1;
    else  % Вторая область взаимодействия
        R_min = R_min2;
        R_max = R_max2;
    end

    % Рассчет взаимодействия через мемристор
    V = abs(y(1) - y(4)); % Разность напряжений двух генераторов
    R_memristor = calculate_R_memristor_multi_region(V, R_min, R_max);
    I_coupling = V / R_memristor; % Ток через мемристор

    % Модель двух связанных генераторов
    dy = zeros(6, 1); % Вектор производных

    % Первый генератор
    dy(1:3) = model_FHN(t, y(1:3), etta1, epsilon1);
    dy(1) = dy(1); % Влияние мемристивного тока

    % Второй генератор
    dy(4:6) = model_FHN(t, y(4:6), etta2, epsilon2);
    dy(4) = dy(4) + I_coupling; % Влияние мемристивного тока
end

function R = calculate_R_memristor_multi_region(V, R_min, R_max)
    % Модель нелинейного сопротивления мемристора
    R = R_min + (R_max - R_min) * exp(abs(V)); % Нелинейное изменение
end

function dy1 = model_FHN(t1, y1, etta1, epsilon1)
    % Модель одного генератора ФицХью-Нагумо

    Fi_ion_1 = 0.8; % eV потенциальный барьер

    q = 1.6e-19;   % Заряд электрона
    k = 1.38e-23;  % Постоянная Больцмана
    T = 300;       % Температура
    S = 5;         % Площадь

    % Параметры
    A_1 = 2.8e3 * exp(1 - 300 / T);
    A_2 = -2.8e3 * exp(1 - 300 / T);
    B_1 = 5.3e21;
    B_2 = 5.3e21;
    Vset_1 = -1;
    Vreset_1 = 2;

    dy1 = zeros(3, 1);

    if y1(1) >= 0
        j_1_lin = A_1 * y1(1);
        j_1_nonlin = B_1 * y1(1) * exp(-44.04 + 1.3 * sqrt(y1(1)));

    else
        j_1_lin = A_2 * abs(y1(1));
        j_1_nonlin = B_2 * abs(y1(1)) * exp(-44.04 + 1.3 * sqrt(abs(y1(1))));
    end

    J_FHN = 0.821 * (y1(3) * j_1_lin + (1 - y1(3)) * j_1_nonlin)*1000 * S * 1e-12;

    dy1(1) = J_FHN - y1(2) + 0.008 * sin(pi * y1(1));
    dy1(2) = epsilon1 * (g(y1(1)) - y1(2) - etta1);

    if y1(1) > Vset_1
        dy1(3) = +1e8 * (exp(q * (-Fi_ion_1 - 0.2 * y1(1)) / (k * T)) * (1 - (2 * y1(3) - 1)^20));
    elseif y1(1) < Vreset_1
        dy1(3) = -1e8 * (exp(q * (-Fi_ion_1 + 0.2 * y1(1)) / (k * T)) * (1 - (2 * y1(3) - 1)^20));
    else
        dy1(3) = 0;
    end

    function G = g(y)
        al_1 = 0.5;
        be_1 = 1.5;
        if y < 0
            G = al_1 * y;
        else
            G = be_1 * y;
        end
    end
end

function period_error = calculate_period_error(t, y)
    % Вычисление ошибки синхронизации через разницу периодов генераторов
    voltage1 = y(:, 1); % Напряжение первого генератора
    voltage2 = y(:, 4); % Напряжение второго генератора

    % Найдём моменты пересечения уровня 0 в положительном направлении
    crossings1 = find_crossings(t, voltage1);
    crossings2 = find_crossings(t, voltage2);

    % Рассчёт периодов
    periods1 = diff(crossings1);
    periods2 = diff(crossings2);

    % Сравниваем периоды генераторов
    min_cycles = min(length(periods1), length(periods2));
    period_error = (periods1(1:min_cycles) - periods2(1:min_cycles));
end

function crossings = find_crossings(t, signal)
    % Поиск пересечений уровня 0 в положительном направлении
    crossings = [];
    for i = 1:length(signal)-1
        if signal(i) < 0 && signal(i+1) >= 0 % Пересечение уровня 0
            crossing_time = t(i) + (t(i+1) - t(i)) * (-signal(i) / (signal(i+1) - signal(i)));
            crossings = [crossings; crossing_time];
        end
    end
end
