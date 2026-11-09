%% Анализ вычислительных экспериментов метода Ричардсона
% Поместите этот файл в одну папку с тремя CSV-файлами и запустите.
% Скрипт строит графики и сохраняет их в папку figures.

clear; clc; close all;

files.eps   = '../data/experiment_eps.csv';
files.det   = '../data/experiment_det.csv';
files.order = '../data/experiment_order.csv';

% Проверка наличия исходных файлов
names = fieldnames(files);
for k = 1:numel(names)
    f = files.(names{k});
    if ~isfile(f)
        error('Не найден файл %s. Поместите CSV рядом со скриптом.', f);
    end
end

outDir = 'figures';
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% Чтение таблиц. VariableNamingRule сохраняет имена столбцов CSV.
E = readtable(files.eps,   'VariableNamingRule', 'preserve');
D = readtable(files.det,   'VariableNamingRule', 'preserve');
O = readtable(files.order, 'VariableNamingRule', 'preserve');

%% Опыт 1. Изменение требуемой точности eps
epsVal = E.("eps");
[epsVal, ix] = sort(epsVal, 'descend');
E = E(ix,:);

figure('Name', 'Опыт 1 — точность', 'Color', 'w');
tiledlayout(2,1);

nexttile;
semilogx(epsVal, E.cycles, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('\epsilon');
ylabel('Число полных циклов');
title('Зависимость числа циклов от требуемой точности');

nexttile;
loglog(epsVal, E.actual_error, 'o-', 'LineWidth', 1.5); hold on;
loglog(epsVal, E.error_estimate, 's--', 'LineWidth', 1.5);
loglog(epsVal, E.residual_norm, 'o-', 'LineWidth', 1.5);
loglog(epsVal, epsVal, 'k:', 'LineWidth', 1.3);
grid on;
xlabel('\epsilon');
ylabel('Ошибка');
legend('Фактическая ошибка ||x-x^*||', ...
       'Оценка ошибки по критерию остановки', ...
       'Невязка', ...
       '\epsilon', ...
       'Location', 'best');
title('Фактическая ошибка и оценка ошибки');
saveas(gcf, fullfile(outDir, '01_epsilon2.png')); 

%% Опыт 2. Изменение определителя матрицы D
detVal = D.("det_D");
[detVal, ix] = sort(detVal, 'descend');
D = D(ix,:);

% Оценка обусловленности диагональной положительной матрицы D
kappaD = D.("Delta") ./ D.("delta");

figure('Name', 'Опыт 2 — определитель', 'Color', 'w');
tiledlayout(2,2);

nexttile;
semilogx(detVal, D.cycles, 'o-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('det(D)');
ylabel('Число полных циклов');
title('Число циклов');

nexttile;
loglog(detVal, D.("delta"), 'o-', 'LineWidth', 1.5); hold on;
loglog(detVal, D.("Delta"), 's-', 'LineWidth', 1.5);
grid on;
xlabel('det(D)');
ylabel('Граница спектра');
legend('\delta', '\Delta', 'Location', 'best');
title('Спектральные границы');

nexttile;
loglog(detVal, kappaD, 'd-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('det(D)');
ylabel('\kappa_2(D) = \Delta/\delta');
title('Число обусловленности D');

nexttile;
semilogx(detVal, D.M3, '^-', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('det(D)');
ylabel('M_3');
title('Коэффициент оценки сходимости');
saveas(gcf, fullfile(outDir, '02_determinant.png'));

figure('Name', 'Опыт 2 — ошибки', 'Color', 'w');
loglog(detVal, D.actual_error, 'o-', 'LineWidth', 1.5); hold on;
loglog(detVal, D.error_estimate, 's--', 'LineWidth', 1.5);
loglog(detVal, D.residual_norm, '^-', 'LineWidth', 1.5);
yline(D.eps(1), 'k:', '\epsilon');
grid on;
xlabel('det(D)');
ylabel('Ошибка');
legend('Фактическая ошибка', 'Оценка ошибки', 'Невязка', '\epsilon', ...
       'Location', 'best');
title('Ошибки решения при изменении det(D)');
saveas(gcf, fullfile(outDir, '03_determinant_errors.png'));

%% Опыт 3. Перестановка трех параметров alpha
% Порядок задается строками вида "0-1-2".
labels = string(O.("order"));
x = categorical(labels, labels, 'Ordinal', true);

figure('Name', 'Опыт 3 — порядок параметров', 'Color', 'w');
tiledlayout(2,1);

nexttile;
bar(x, O.cycles);
grid on;
ylabel('Число полных циклов');
title('Число циклов для разных порядков параметров');
xlabel('Порядок индексов \alpha');

nexttile;
bar(x, [O.first_cycle_error_step1, ...
        O.first_cycle_error_step2, ...
        O.first_cycle_error_step3]);
grid on;
ylabel('||x^{(k)}-x^*||');
xlabel('Порядок индексов \alpha');
legend('После шага 1', 'После шага 2', 'После шага 3', ...
       'Location', 'best');
title('Ошибка на промежуточных шагах первого цикла');
saveas(gcf, fullfile(outDir, '04_parameter_order.png'));

%% Краткая численная сводка в командном окне
fprintf('\nОПЫТ 1 — eps\n');
fprintf('Фактическая точность достигнута в %d из %d запусков.\n', ...
    sum(E.actual_accuracy_met == 1), height(E));
fprintf('Минимальная eps: %.1e, фактическая ошибка: %.3e.\n', ...
    min(E.("eps")), E.actual_error(E.("eps") == min(E.("eps"))));

fprintf('\nОПЫТ 2 — det(D)\n');
fprintf('Число циклов: от %d до %d.\n', min(D.cycles), max(D.cycles));
fprintf('Число обусловленности D: от %.3e до %.3e.\n', ...
    min(kappaD), max(kappaD));

fprintf('\nОПЫТ 3 — порядок alpha\n');
fprintf('Число циклов по всем перестановкам: ');
fprintf('%d ', O.cycles);
fprintf('\n');

fprintf('\nГрафики сохранены в папке "%s".\n', outDir);
