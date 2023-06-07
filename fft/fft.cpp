#include <iostream>
#include <complex>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h> 

class FFT {
public:
    // Конструктор
    FFT() {}

    // Прямое преобразование Фурье
    std::vector<std::complex<double>> fft(std::vector<std::complex<double>>& x) {
        int n = x.size();

        // Базовый полином
        if (n == 1) return x;

        // Распаковка значений из массива
        std::vector<std::complex<double>> even(n / 2), odd(n / 2);
        for (int i = 0; i < n / 2; ++i) {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }

        // Рекурсивные вызовы
        std::vector<std::complex<double>> q = fft(even);
        std::vector<std::complex<double>> r = fft(odd);

        // Сложение результатов
        std::vector<std::complex<double>> y(n);
        for (int i = 0; i < n / 2; ++i) {
            std::complex<double> t = std::polar(1.0, -2 * M_PI * i / n) * r[i];
            y[i] = q[i] + t;
            y[i + n / 2] = q[i] - t;
        }
        return y;
    }

    // Обратное преобразование Фурье
    std::vector<std::complex<double>> ifft(std::vector<std::complex<double>>& x) {
        // Конвертация вектора комплексных чисел к обратному порядку значений
        std::vector<std::complex<double>> y = fft(x);
        std::vector<std::complex<double>> z(y.size());
        for (int i = 0; i < y.size(); ++i) {
            z[i] = y[(y.size() - i) % y.size()] / static_cast<double>(y.size());
        }
        return z;
    }
};

// Пример использования класса
int main() {
    FFT fft;

    // Создание случайного вектора комплексных чисел длиной кратной 2, 3 или 5
    const int N = 20;
    std::vector<std::complex<double>> x(N);
    for (int i = 0; i < x.size(); ++i) {
        double re = (double)rand() / RAND_MAX;
        double im = (double)rand() / RAND_MAX;
        x[i] = std::complex<double>(re, im);
    }

    // Применение прямого и обратного преобразования
    std::vector<std::complex<double>> y = fft.fft(x);
    std::vector<std::complex<double>> z = fft.ifft(y);

    // Вывод результатов
    std::cout << "Input: ";
    for (auto& elem : x) {
        std::cout << elem << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "FFT: ";
    for (auto& elem : y) {
        std::cout << elem << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "Inverse FFT: ";
    for (auto& elem : z) {
        std::cout << elem << " ";
    }
    std::cout << std::endl << std::endl;

    // Вычисление ошибки между входными и выходными данными

    double error_lst[N];
    double error = 0.0;
    std::cout << "Errors: " << std::endl;
    for (int i = 0; i < x.size(); ++i) {
        error_lst[i] = std::abs(x[i] - z[i]);
        error += std::abs(x[i] - z[i]);
        std::cout << error_lst[i] << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "Error: " << error << std::endl;

    return 0;
}