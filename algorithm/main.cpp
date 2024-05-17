#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <vector>

#define M_PI 3.14159265358979323846

// Границы прямоугольной области
const double a = 0.0;
const double b = 1.0;
const double c = 0.0;
const double d = 1.0;

// Вариант - 1
// Основная задача
double f_main(double x, double y)
{
    double temp = sin(M_PI * x * y);
    return temp * temp;
}

// Основные граничные условия
double mu1_main(double y)
{
    return sin(M_PI * y);
}
double mu2_main(double y)
{
    return sin(M_PI * y);
}
double mu3_main(double x)
{
    return x - x * x;
}
double mu4_main(double x)
{
    return x - x * x;
}

// Функция u для тестовой задачи
double u_test(double x, double y)
{
    return exp(f_main(x, y));
}
// Функция f для тестовой задачи
double f_test(double x, double y)
{
    return 0.5 * u_test(x, y) * M_PI * M_PI * (x * x + y * y) * (-1 - 4 * cos(2 * M_PI * x * y) + cos(4 * M_PI * x * y));
}

// Подставляем 0 и 1
// Тестовые граничные условия
double mu1_test(double y)
{
    return 1;
}
double mu2_test(double y)
{
    double temp = sin(M_PI * y);
    return exp(temp * temp);
}
double mu3_test(double x)
{
    return 1;
}
double mu4_test(double x)
{
    double temp = sin(M_PI * x);
    return exp(temp * temp);
}

// Совокупность граничных функций
double gran_func_test(double x, double y, int i, int j, int n, int m)
{
    if (i == 0)
        return mu1_test(y);
    if (i == n)
        return mu2_test(y);
    if (j == 0)
        return mu3_test(x);
    if (j == m)
        return mu4_test(x);
}
double gran_func_main(double x, double y, int i, int j, int n, int m)
{
    if (i == 0)
        return mu1_main(y);
    if (i == n)
        return mu2_main(y);
    if (j == 0)
        return mu3_main(x);

    if (j == m)
        return mu4_main(x);
}

// Решение тестовой задачи
extern "C" __declspec(dllexport) void solve_test(int n, int m, int n_max, double eps)
{
    std::ofstream outputFile_numerical("result//numerical.txt");
    std::ofstream outputFile_auxiliary("result//auxiliary.txt");
    std::ofstream outputFile_diff_solutions("result//diff_solutions.txt");
    std::ofstream outputFile_results("result//results.txt");

    // Задаём шаг по x и y
    double h = (b - a) / n;
    double k = (d - c) / m;

    // Создаём матрицу численного решения и матрицу невязки
    std::vector<std::vector<double>> v, r_s;

    // Результат выхода из программы
    int result;

    // Количество итераций
    int S = 0;

    // Задаём максимальную погрешность и коэффициенты
    double epsMax = 0.0;
    double diff = 0.0;
    double h2 = ((n / (b - a)) * (n / (b - a)));
    double k2 = ((m / (d - c)) * (m / (d - c)));
    double a2 = -2 * (h2 + k2);

    // Создаём вектора сеток по x и y
    std::vector<double> x(n + 1);
    std::vector<double> y(m + 1);

    // Заполняем сетки
    x[0] = a;
    for (int i = 1; i <= n; ++i)
    {
        x[i] = a + i * h;
    }

    if (outputFile_numerical.is_open())
    {
        for (int i = 0; i < n + 1; ++i)
        {
            outputFile_numerical << x[i] << ' ';
        }
        outputFile_numerical << "\n\n";
    }

    y[0] = c;
    for (int j = 1; j <= m; ++j)
    {
        y[j] = c + j * k;
    }

    if (outputFile_numerical.is_open())
    {
        for (int j = 0; j < m + 1; ++j)
        {
            outputFile_numerical << y[j] << ' ';
        }
        outputFile_numerical << "\n\n";
    }

    // Создаём промежуточный вектор для заполнения
    std::vector<double> vec(m + 1, 0);
    for (int t = 0; t <= m; t++)
    {
        vec.push_back(0);
    }

    // Заполняем численное решение и невязку
    for (int i = 0; i <= n; ++i)
    {
        v.push_back(vec);
        r_s.push_back(vec);
    }

    // Подставляем граничные условия
    for (int j = 0; j <= m; j++)
    {
        for (int i = 0; i <= n; i++)
        {
            if ((i == 0) || (i == n) || (j == 0) || (j == m))
            {
                v[i][j] = gran_func_test(x[i], y[j], i, j, n, m);
            }
        }
    }

    // Норма невязки
    double current_norm = 0.0;
    // Параметр, по которому определяется численное решение
    double tau = 0.0;

    // Для замера времени (начало)
    clock_t start = clock();

    while (true)
    {
        if (S >= n_max)
        {
            if (outputFile_numerical.is_open() && outputFile_diff_solutions.is_open() && outputFile_auxiliary.is_open())
            {
                for (int j = 0; j < m + 1; ++j)
                {
                    for (int i = 0; i < n + 1; ++i)
                    {
                        outputFile_numerical << v[i][j] << ' ';
                        outputFile_auxiliary << u_test(x[i], y[j]) << ' ';
                        outputFile_diff_solutions << u_test(x[i], y[j]) - v[i][j] << ' ';
                    }
                    outputFile_numerical << '\n';
                    outputFile_auxiliary << '\n';
                    outputFile_diff_solutions << '\n';
                }
                outputFile_numerical << '\n';
                outputFile_auxiliary << '\n';
                outputFile_diff_solutions << '\n';
            }

            result = 0;

            if (outputFile_results.is_open())
            {
                outputFile_results << result << '\n';
                outputFile_results << epsMax << '\n';
                outputFile_results << diff << '\n';
                outputFile_results << S << '\n';
            }

            break;
        }

        // Обновление переменных
        epsMax = 0.0;
        diff = 0.0;
        current_norm = 0.0;

        // Высчитываем невязку и её норму по бесконечности
        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < m; j++)
            {
                r_s[i][j] = -a2 * v[i][j] - h2 * (v[i + 1][j] + v[i - 1][j]) - k2 * (v[i][j + 1] + v[i][j - 1]) - f_test(x[i], y[j]);
                current_norm = fmax(current_norm, fabs(r_s[i][j]));
            }
        }

        // Вспомогательная переменная и переменные для числителя и знаменателя
        double temp = 0.0;
        double s_up = 0.0;
        double s_down = 0.0;

        // Вычисление коэффициента tau
        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < m; j++)
            {
                temp = (r_s[i][j] * -a2 - h2 * (r_s[i + 1][j] + r_s[i - 1][j]) - k2 * (r_s[i][j + 1] + r_s[i][j - 1]));
                s_up += temp * r_s[i][j];
                s_down += temp * temp;
            }
        }
        tau = s_up / s_down;

        // Вспомогательная переменная для погрешности
        double local_epsMax = epsMax;
        double local_diff = diff;

        // Обновление численного решения и вычисление погрешности
        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < m; j++)
            {
                double v_old_local = v[i][j];
                double v_new_local = v[i][j] - tau * r_s[i][j];
                local_epsMax = fmax(local_epsMax, fabs(v_old_local - v_new_local));
                v[i][j] = v_new_local;
                local_diff = fmax(local_diff, fabs(u_test(x[i], y[j]) - v_new_local));
            }
        }
        epsMax = local_epsMax;
        diff = local_diff;

        ++S;

        if (epsMax < eps)
        {
            if (outputFile_numerical.is_open() && outputFile_diff_solutions.is_open())
            {
                for (int j = 0; j < m + 1; ++j)
                {
                    for (int i = 0; i < n + 1; ++i)
                    {
                        outputFile_numerical << v[i][j] << ' ';
                        outputFile_auxiliary << u_test(x[i], y[j]) << ' ';
                        outputFile_diff_solutions << u_test(x[i], y[j]) - v[i][j] << ' ';
                    }
                    outputFile_numerical << '\n';
                    outputFile_auxiliary << '\n';
                    outputFile_diff_solutions << '\n';
                }
                outputFile_numerical << '\n';
                outputFile_auxiliary << '\n';
                outputFile_diff_solutions << '\n';
            }

            result = 1;

            if (outputFile_results.is_open())
            {
                outputFile_results << result << '\n';
                outputFile_results << epsMax << '\n';
                outputFile_results << diff << '\n';
                outputFile_results << S << '\n';
            }

            break;
        }
    }

    // Для замера времени (конец)
    clock_t stop = clock();
    double seconds = (double)(stop - start) / CLOCKS_PER_SEC;

    if (outputFile_results.is_open())
    {
        outputFile_results << seconds << "\n\n";
    }

    outputFile_numerical.close();
    outputFile_auxiliary.close();
    outputFile_diff_solutions.close();
    outputFile_results.close();
}

extern "C" __declspec(dllexport) void solve_main(int n, int m, int n_max, double eps, int n_max_, double eps_)
{
    std::ofstream outputFile_numerical("result//numerical.txt");
    std::ofstream outputFile_diff_solutions("result//diff_solutions.txt");
    std::ofstream outputFile_results("result//results.txt");

    int n_ = 2 * n;
    int m_ = 2 * m;

    double h = (b - a) / n;
    double k = (d - c) / m;
    double h_ = (b - a) / n_;
    double k_ = (d - c) / m_;

    std::vector<std::vector<double>> v, r_s;
    std::vector<std::vector<double>> v_, r_s_;

    int result;
    int result_;

    int S = 0;
    int S_ = 0;

    double precision = 0.0;

    double epsMax = 0;
    double epsMax_ = 0;

    double h2 = ((n / (b - a)) * (n / (b - a)));
    double k2 = ((m / (d - c)) * (m / (d - c)));
    double h2_ = ((n_ / (b - a)) * (n_ / (b - a)));
    double k2_ = ((m_ / (d - c)) * (m_ / (d - c)));

    double a2 = -2 * (h2 + k2);
    double a2_ = -2 * (h2_ + k2_);

    std::vector<double> x(n + 1);
    std::vector<double> y(m + 1);
    std::vector<double> x_(n_ + 1);
    std::vector<double> y_(m_ + 1);

    x[0] = a;
    for (int i = 1; i <= n; ++i)
    {
        x[i] = a + i * h;
    }

    if (outputFile_numerical.is_open())
    {
        for (int i = 0; i < n + 1; ++i)
        {
            outputFile_numerical << x[i] << ' ';
        }
        outputFile_numerical << "\n\n";
    }

    y[0] = c;
    for (int j = 1; j <= m; ++j)
    {
        y[j] = c + j * k;
    }

    if (outputFile_numerical.is_open())
    {
        for (int j = 0; j < m + 1; ++j)
        {
            outputFile_numerical << y[j] << ' ';
        }
        outputFile_numerical << "\n\n";
    }

    x_[0] = a;
    for (int i = 1; i <= n_; ++i)
    {
        x_[i] = a + i * h_;
    }

    if (outputFile_numerical.is_open())
    {
        for (int i = 0; i < n_ + 1; ++i)
        {
            outputFile_numerical << x_[i] << ' ';
        }
        outputFile_numerical << "\n\n";
    }

    y_[0] = c;
    for (int j = 1; j <= m_; ++j)
    {
        y_[j] = c + j * k_;
    }

    if (outputFile_numerical.is_open())
    {
        for (int j = 0; j < m_ + 1; ++j)
        {
            outputFile_numerical << y_[j] << ' ';
        }
        outputFile_numerical << "\n\n";
    }

    std::vector<double> vec(m + 1, 0);
    for (int t = 0; t <= m; t++)
    {
        vec.push_back(0);
    }
    std::vector<double> vec_(m_ + 1, 0);
    for (int t = 0; t <= m_; t++)
    {
        vec_.push_back(0);
    }

    for (int i = 0; i <= n; ++i)
    {
        v.push_back(vec);
        r_s.push_back(vec);
    }
    for (int i = 0; i <= n_; ++i)
    {
        v_.push_back(vec_);
        r_s_.push_back(vec_);
    }

    for (int j = 0; j <= m; j++)
    {
        for (int i = 0; i <= n; i++)
        {
            if ((i == 0) || (i == n) || (j == 0) || (j == m))
            {
                v[i][j] = gran_func_main(x[i], y[j], i, j, n, m);
            }
        }
    }
    for (int j = 0; j <= m_; j++)
    {
        for (int i = 0; i <= n_; i++)
        {
            if ((i == 0) || (i == n_) || (j == 0) || (j == m_))
            {
                v_[i][j] = gran_func_main(x_[i], y_[j], i, j, n_, m_);
            }
        }
    }

    double current_norm = 0.0;
    double current_norm_ = 0.0;
    double tau = 0.0;
    double tau_ = 0.0;

    clock_t start = clock();
    clock_t stop;

    while (true)
    {
        if (S >= n_max)
        {
            if (outputFile_numerical.is_open())
            {
                for (int j = 0; j < m + 1; ++j)
                {
                    for (int i = 0; i < n + 1; ++i)
                    {
                        outputFile_numerical << v[i][j] << ' ';
                    }
                    outputFile_numerical << '\n';
                }
                outputFile_numerical << '\n';
            }

            result = 0;
            stop = clock();
            double seconds = (double)(stop - start) / CLOCKS_PER_SEC;

            if (outputFile_results.is_open())
            {
                outputFile_results << result << '\n';
                outputFile_results << seconds << '\n';
                outputFile_results << S << '\n';
                outputFile_results << epsMax << "\n\n";
            }

            break;
        }

        epsMax = 0.0;
        current_norm = 0.0;

        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < m; j++)
            {
                r_s[i][j] = -a2 * v[i][j] - h2 * (v[i + 1][j] + v[i - 1][j]) - k2 * (v[i][j + 1] + v[i][j - 1]) - f_main(x[i], y[j]);
                current_norm = fmax(current_norm, fabs(r_s[i][j]));
            }
        }

        double temp = 0.0;
        double s_up = 0.0;
        double s_down = 0.0;

        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < m; j++)
            {
                temp = (r_s[i][j] * -a2 - h2 * (r_s[i + 1][j] + r_s[i - 1][j]) - k2 * (r_s[i][j + 1] + r_s[i][j - 1]));
                s_up += temp * r_s[i][j];
                s_down += temp * temp;
            }
        }
        tau = s_up / s_down;

        double local_epsMax = epsMax;

        for (int i = 1; i < n; i++)
        {
            for (int j = 1; j < m; j++)
            {
                double v_old_local = v[i][j];
                double v_new_local = v[i][j] - tau * r_s[i][j];
                v[i][j] = v_new_local;
                local_epsMax = fmax(local_epsMax, fabs(v_old_local - v_new_local));
            }
        }
        epsMax = local_epsMax;

        ++S;

        if (epsMax < eps)
        {
            if (outputFile_numerical.is_open())
            {
                for (int j = 0; j < m + 1; ++j)
                {
                    for (int i = 0; i < n + 1; ++i)
                    {
                        outputFile_numerical << v[i][j] << ' ';
                    }
                    outputFile_numerical << '\n';
                }
                outputFile_numerical << '\n';
            }

            result = 1;
            stop = clock();
            double seconds = (double)(stop - start) / CLOCKS_PER_SEC;

            if (outputFile_results.is_open())
            {
                outputFile_results << result << '\n';
                outputFile_results << seconds << '\n';
                outputFile_results << S << '\n';
                outputFile_results << epsMax << "\n\n";
            }

            break;
        }
    }

    clock_t start_ = clock();

    while (true)
    {
        if (S_ >= n_max_)
        {
            if (outputFile_numerical.is_open())
            {
                for (int j = 0; j < m_ + 1; ++j)
                {
                    for (int i = 0; i < n_ + 1; ++i)
                    {
                        outputFile_numerical << v_[i][j] << ' ';
                    }
                    outputFile_numerical << '\n';
                }
                outputFile_numerical << "\n\n";
            }

            result_ = 0;

            if (outputFile_results.is_open())
            {
                outputFile_results << result_ << '\n';
                outputFile_results << S_ << '\n';
                outputFile_results << epsMax_ << '\n';
            }

            break;
        }

        epsMax_ = 0.0;
        current_norm_ = 0.0;

        for (int i = 1; i < n_; i++)
        {
            for (int j = 1; j < m_; j++)
            {
                r_s_[i][j] = -a2_ * v_[i][j] - h2_ * (v_[i + 1][j] + v_[i - 1][j]) - k2_ * (v_[i][j + 1] + v_[i][j - 1]) - f_main(x_[i], y_[j]);
                current_norm_ = fmax(current_norm_, fabs(r_s_[i][j]));
            }
        }

        double temp_ = 0.0;
        double s_up_ = 0.0;
        double s_down_ = 0.0;

        for (int i = 1; i < n_; i++)
        {
            for (int j = 1; j < m_; j++)
            {
                temp_ = (r_s_[i][j] * -a2_ - h2_ * (r_s_[i + 1][j] + r_s_[i - 1][j]) - k2_ * (r_s_[i][j + 1] + r_s_[i][j - 1]));
                s_up_ += temp_ * r_s_[i][j];
                s_down_ += temp_ * temp_;
            }
        }
        tau_ = s_up_ / s_down_;

        double local_epsMax_ = epsMax_;

        for (int i = 1; i < n_; i++)
        {
            for (int j = 1; j < m_; j++)
            {
                double v_old_local_ = v_[i][j];
                double v_new_local_ = v_[i][j] - tau_ * r_s_[i][j];
                v_[i][j] = v_new_local_;
                local_epsMax_ = fmax(local_epsMax_, fabs(v_old_local_ - v_new_local_));
            }
        }
        epsMax_ = local_epsMax_;

        ++S_;

        if (epsMax_ < eps_)
        {
            if (outputFile_numerical.is_open())
            {
                for (int j = 0; j < m_ + 1; ++j)
                {
                    for (int i = 0; i < n_ + 1; ++i)
                    {
                        outputFile_numerical << v_[i][j] << ' ';
                    }
                    outputFile_numerical << '\n';
                }
                outputFile_numerical << "\n\n";
            }

            result_ = 1;

            if (outputFile_results.is_open())
            {
                outputFile_results << result_ << '\n';
                outputFile_results << S_ << '\n';
                outputFile_results << epsMax_ << '\n';
            }

            break;
        }
    }

    clock_t stop_ = clock();

    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            precision = fmax(precision, abs(v[i][j] - v_[2 * i][2 * j]));
        }
    }

    if (outputFile_diff_solutions.is_open())
    {
        for (int j = 0; j < m + 1; j++)
        {
            for (int i = 0; i < n + 1; i++)
            {
                outputFile_diff_solutions << v[i][j] - v_[2 * i][2 * j] << ' ';
            }
            outputFile_diff_solutions << '\n';
        }
        outputFile_diff_solutions << '\n';
    }

    double seconds_ = (double)(stop_ - start_) / CLOCKS_PER_SEC;

    if (outputFile_results.is_open())
    {
        outputFile_results << precision << '\n';
        outputFile_results << seconds_ << "\n\n";
    }

    outputFile_numerical.close();
    outputFile_diff_solutions.close();
    outputFile_results.close();
}