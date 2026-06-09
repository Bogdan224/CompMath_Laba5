using System;
using System.Collections.Generic;
using System.Linq;
using ScottPlot;

namespace Lab5Smoothing
{
    class Program
    {
        static Random rand = new Random();

        static void Main(string[] args)
        {
            // 1. Генерация исходных данных
            int n = 200; // количество точек от 0 до 199
            double[] x = new double[n];
            double[] y_raw = new double[n];

            Console.WriteLine("Первые 10 исходных точек");
            Console.WriteLine("i\tsqrt(i/10)\trnd(1)\tYi");
            for (int i = 0; i < n; i++)
            {
                x[i] = i;
                double trend = Math.Sqrt(i / 10.0);
                double noise = rand.NextDouble(); // rnd(1) — [0,1)
                y_raw[i] = trend + noise;
                if(i < 10)
                    Console.WriteLine($"{i}\t{trend:F3}\t\t{noise:F4}\t{y_raw[i]:F4}");
            }

            // 2. Методы сглаживания
            double[] y_sma3 = MovingAverage(y_raw, 3);   // скользящее среднее окно 3
            double[] y_sma7 = MovingAverage(y_raw, 7);   // скользящее среднее окно 7
            double[] y_localLin = LocalPoly1(x, y_raw, 3); // локальный МНК (степень 1, окно 3)
            double[] y_globalQuad = GlobalPoly2(x, y_raw); // глобальный МНК (степень 2)

            Plot plt = new Plot();

            // Исходные данные (серые точки)
            var scatterRaw = plt.Add.Scatter(x, y_raw);
            scatterRaw.Color = Colors.Gray;
            scatterRaw.MarkerSize = 3;
            scatterRaw.LineWidth = 0;
            scatterRaw.Label = "Исходные данные (зашумленные)";

            // Скользящее среднее окно 3 (синий)
            var scatterSma3 = plt.Add.Scatter(x, y_sma3);
            scatterSma3.Color = Colors.Blue;
            scatterSma3.LineWidth = 2;
            scatterSma3.MarkerSize = 0;
            scatterSma3.Label = "Скользящее среднее (окно 3)";

            // Скользящее среднее окно 7 (красный)
            var scatterSma7 = plt.Add.Scatter(x, y_sma7);
            scatterSma7.Color = Colors.Red;
            scatterSma7.LineWidth = 2;
            scatterSma7.MarkerSize = 0;
            scatterSma7.Label = "Скользящее среднее (окно 7)";

            // Локальный МНК (зелёный)
            var scatterLocal = plt.Add.Scatter(x, y_localLin);
            scatterLocal.Color = Colors.Green;
            scatterLocal.LineWidth = 2;
            scatterLocal.MarkerSize = 0;
            scatterLocal.Label = "Локальный МНК (степень 1, окно 3)";

            // Глобальный МНК (чёрный)
            var scatterGlobal = plt.Add.Scatter(x, y_globalQuad);
            scatterGlobal.Color = Colors.Black;
            scatterGlobal.LineWidth = 2;
            scatterGlobal.MarkerSize = 0;
            scatterGlobal.Label = "Глобальный МНК (степень 2)";

            // Настройка графика
            plt.Title("Сравнение методов сглаживания");
            plt.XLabel("Индекс i (Xi = i)");
            plt.YLabel("Значение Yi");
            plt.ShowLegend(Alignment.UpperLeft);

            // Добавление сетки
            plt.Grid.MajorLineColor = Colors.LightGray;
            plt.Grid.MajorLinePattern = LinePattern.Solid;

            // Сохранение графика
            string path = $"C:\\Users\\{Environment.UserName}\\Downloads\\Laba5\\";
            if(!Directory.Exists(path))
                Directory.CreateDirectory(path);
            string filename = "smoothing_results.png";
            plt.SavePng(path + filename, 800, 600);
            Console.WriteLine($"\nГрафик сохранён в файл: {path}{filename}");

            // 4. Вывод параметров тренда
            Console.WriteLine("\nРезультаты глобального МНК (полином 2 степени)");
            (double a0, double a1, double a2) = FitPoly2(x, y_raw);
            Console.WriteLine($"P2(x) = {a0:F6} + {a1:F6}*x + {a2:F10}*x^2");
            Console.WriteLine($"Теоретический тренд: sqrt(x/10)");

            Console.WriteLine("\nПрограмма завершена. Нажмите любую клавишу для выхода...");
            Console.ReadKey();
        }

        /// <summary>
        /// Метод скользящего среднего
        /// </summary>
        static double[] MovingAverage(double[] data, int windowSize)
        {
            if (windowSize % 2 == 0)
                throw new ArgumentException("Размер окна должен быть нечётным");

            int n = data.Length;
            double[] result = new double[n];
            int half = windowSize / 2;

            for (int i = 0; i < n; i++)
            {
                int start = Math.Max(0, i - half);
                int end = Math.Min(n - 1, i + half);
                double sum = 0;
                int count = 0;

                for (int j = start; j <= end; j++)
                {
                    sum += data[j];
                    count++;
                }

                result[i] = sum / count;
            }

            return result;
        }

        /// <summary>
        /// Локальный МНК с полиномом 1 степени (линейная аппроксимация) по окну соседних точек
        /// </summary>
        static double[] LocalPoly1(double[] x, double[] y, int windowSize)
        {
            if (windowSize % 2 == 0)
                throw new ArgumentException("Размер окна должен быть нечётным");

            int n = x.Length;
            double[] result = new double[n];
            int half = windowSize / 2;

            for (int i = 0; i < n; i++)
            {
                int start = Math.Max(0, i - half);
                int end = Math.Min(n - 1, i + half);
                int m = end - start + 1; // фактический размер окна

                // Подготовка данных для локального МНК
                double[] xLocal = new double[m];
                double[] yLocal = new double[m];
                for (int j = 0; j < m; j++)
                {
                    xLocal[j] = x[start + j];
                    yLocal[j] = y[start + j];
                }

                // Решаем СЛАУ для полинома 1 степени: a0 + a1*x
                // Нормальные уравнения:
                // [ m     sum(x) ] [a0] = [ sum(y)    ]
                // [ sum(x) sum(x^2)] [a1]   [ sum(x*y) ]
                double sumX = xLocal.Sum();
                double sumY = yLocal.Sum();
                double sumX2 = xLocal.Select(v => v * v).Sum();
                double sumXY = xLocal.Zip(yLocal, (xv, yv) => xv * yv).Sum();

                // Решение методом Гаусса для 2 переменных
                double[,] A = { { m, sumX }, { sumX, sumX2 } };
                double[] B = { sumY, sumXY };
                double[] coeff = SolveGauss(A, B);

                // Вычисляем значение в центральной точке
                result[i] = coeff[0] + coeff[1] * x[i];
            }

            return result;
        }

        /// <summary>
        /// Глобальный МНК с полиномом 2 степени (парабола) по всем точкам
        /// </summary>
        static double[] GlobalPoly2(double[] x, double[] y)
        {
            var coeff = FitPoly2(x, y);
            double a0 = coeff.a0, a1 = coeff.a1, a2 = coeff.a2;

            double[] result = new double[x.Length];
            for (int i = 0; i < x.Length; i++)
            {
                result[i] = a0 + a1 * x[i] + a2 * x[i] * x[i];
            }
            return result;
        }

        /// <summary>
        /// Подгонка полинома 2 степени (a0 + a1*x + a2*x^2) методом наименьших квадратов
        /// </summary>
        static (double a0, double a1, double a2) FitPoly2(double[] x, double[] y)
        {
            int n = x.Length;
            double sumX = x.Sum();
            double sumX2 = x.Select(v => v * v).Sum();
            double sumX3 = x.Select(v => v * v * v).Sum();
            double sumX4 = x.Select(v => v * v * v * v).Sum();
            double sumY = y.Sum();
            double sumXY = x.Zip(y, (xv, yv) => xv * yv).Sum();
            double sumX2Y = x.Zip(y, (xv, yv) => xv * xv * yv).Sum();

            // Система нормальных уравнений для параболы
            // | n    sumX   sumX2 | |a0|   | sumY   |
            // | sumX sumX2  sumX3 | |a1| = | sumXY  |
            // | sumX2 sumX3 sumX4 | |a2|   | sumX2Y |
            double[,] A = {
                { n, sumX, sumX2 },
                { sumX, sumX2, sumX3 },
                { sumX2, sumX3, sumX4 }
            };
            double[] B = { sumY, sumXY, sumX2Y };

            double[] coeff = SolveGauss(A, B);
            return (coeff[0], coeff[1], coeff[2]);
        }

        /// <summary>
        /// Решение СЛАУ методом Гаусса с выбором главного элемента
        /// </summary>
        static double[] SolveGauss(double[,] A, double[] B)
        {
            int n = B.Length;
            double[,] augmented = new double[n, n + 1];

            // Формируем расширенную матрицу
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                    augmented[i, j] = A[i, j];
                augmented[i, n] = B[i];
            }

            // Прямой ход с выбором главного элемента
            for (int k = 0; k < n; k++)
            {
                // Поиск максимального элемента в столбце k
                int maxRow = k;
                double maxVal = Math.Abs(augmented[k, k]);
                for (int i = k + 1; i < n; i++)
                {
                    if (Math.Abs(augmented[i, k]) > maxVal)
                    {
                        maxVal = Math.Abs(augmented[i, k]);
                        maxRow = i;
                    }
                }

                // Перестановка строк
                if (maxRow != k)
                {
                    for (int j = k; j <= n; j++)
                    {
                        double temp = augmented[k, j];
                        augmented[k, j] = augmented[maxRow, j];
                        augmented[maxRow, j] = temp;
                    }
                }

                if (Math.Abs(augmented[k, k]) < 1e-10)
                    throw new Exception();

                // Нормализация строки k
                for (int j = k + 1; j <= n; j++)
                    augmented[k, j] /= augmented[k, k];
                augmented[k, k] = 1;

                // Вычитание из нижних строк
                for (int i = k + 1; i < n; i++)
                {
                    double factor = augmented[i, k];
                    for (int j = k; j <= n; j++)
                        augmented[i, j] -= factor * augmented[k, j];
                }
            }

            // Обратный ход
            double[] X = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                X[i] = augmented[i, n];
                for (int j = i + 1; j < n; j++)
                    X[i] -= augmented[i, j] * X[j];
            }

            return X;
        }
    }
}