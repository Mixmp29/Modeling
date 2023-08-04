#include <math.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#define _USE_MATH_DEFINES

double generateRandomNumber() {
  static std::random_device rd;
  static std::mt19937 mt(rd());
  static std::uniform_real_distribution<double> dist(0, 1);
  return dist(mt);
}

double kCoef(int n) {
  double result = 1 + 3.322 * log10(n);
  if (result < 5) {
    result = 5;
  }
  if (result > 15) {
    result = 15;
  }
  return result;
}

inline double hStep(double xMax, double xMin, int k) {
  return (xMax - xMin) / k;
}

/* inline double uiCalc(double xi, double mathExpectation,
                     double standartDeviation) {
  return (xi - mathExpectation) / standartDeviation;
}
inline double fiCalc(double u) {
  return pow(M_E, -pow(u, 2) / 2) / sqrt(2 * M_PI);
} */

void calcXiNi(double xMin, double xMax, double n, int k, double h,
              const std::vector<double>& randomNumbers, std::vector<double>& xi,
              std::vector<double>& ni) {
  double left = xMin, right = xMin + h;
  int intervalNumber = 0;

  for (const auto& number : randomNumbers) {
    if (number >= right) {
      xi[intervalNumber] = (left + right) / 2;
      left = right;
      right = left + h;
      ++intervalNumber;
    }
    ++ni[intervalNumber];
  }
  xi[intervalNumber] = (left + right) / 2;
}

void calcWi(int n, int k, const std::vector<double>& ni,
            std::vector<double>& wi) {
  for (int i = 0; i < k; ++i) {
    wi[i] = ni[i] / n;
  }
}

double avrX(int k, std::vector<double>& xi) {
  double avrX = 0;
  for (int i = 0; i < k; i++) {
    avrX += xi[i];
  }

  return avrX / k;
}

double sampleVariance(int n, std::vector<double>& xi, double avrX) {
  double result = 0;
  for (auto& number : xi) {
    result += pow(number - avrX, 2);
  }
  return result / n;
}

double at(int n, std::vector<double>& xi, double avrX, double t) {
  double result = 0;
  double variance = sampleVariance(n, xi, avrX);
  for (int i = 0; i < n - t; i++) {
    result += ((xi[i] - avrX) * (xi[i + t] - avrX)) / (variance * (n - t));
  }

  return result;
}

double calcMathExpectation(int k, std::vector<double>& xi,
                           std::vector<double>& wi) {
  double result = 0;

  for (int i = 0; i < k; i++) {
    result += xi[i] * wi[i];
  }

  return result;
}

double calcDispersion(int k, std::vector<double>& xi, std::vector<double>& wi,
                      double mathExpectation) {
  double result = 0;

  for (int i = 0; i < k; i++) {
    result += pow(xi[i] - mathExpectation, 2) * wi[i];
  }

  return result;
}

double Hi2(int n, int k, double h, std::vector<double>& ni,
           std::vector<double>& ni1, std::vector<double>& xi,
           double mathExpectation, double standartDeviation) {
  double result = 0;
  double a = mathExpectation - std::sqrt(3) * standartDeviation;
  double b = mathExpectation + std::sqrt(3) * standartDeviation;

  std::vector<double> Pi(k);

  // std::cout << "a = " << a << std::endl;
  // std::cout << "b = " << b << std::endl;

  for (int i = 0; i < k - 1; ++i) {
    Pi[i] = xi[i + 1] / (b - a) - xi[i] / (b - a);
    // std::cout << "P[" << i << "] = " << Pi[i] << std::endl;
  }
  Pi[k - 1] = Pi[k - 2];
  // std::cout << "P[" << k - 1 << "] = " << Pi[k - 1] << std::endl;

  for (int i = 0; i < k; ++i) {
    ni1[i] = n * Pi[i];
    // std::cout << "n1[" << i << "] = " << ni1[i] << std::endl;
  }

  for (int i = 0; i < k; ++i) {
    result += std::pow((ni[i] - ni1[i]), 2) / ni1[i];
    // std::cout << "result " << i << " = " << result << std::endl;
  }

  return result;
}

void printInfo(int n, int k, double h, std::vector<double>& ni,
               std::vector<double>& ni1, std::vector<double>& xi,
               std::vector<double>& wi, double mathExpectation,
               double dispersion, double standartDeviation,
               double criterionPirson, double freedomDegree, double averageX) {
  double max = 0;
  int max_i;
  std::ofstream out("Hi2.txt");
  std::ofstream out2("autocorrelation.txt");

  std::cout << n << " случайных чисел в диапазоне от 0 до 1:"
            << "\n";

  std::cout << "+============================+==================+============="
               "======+==="
               "===============+\n";
  std::cout << "|          ai-1; ai          |         ni       |         "
               "xi        |        "
               " wi       |\n";
  std::cout << "+============================+==================+============="
               "======+========"
               "==========+\n";
  for (int i = 0; i < k; ++i) {
    printf(
        "|     %.6f; %.6f     |         %.0f\t|      %lf     |      %lf "
        "   "
        "|\n",
        i * h, i * h + h, ni[i], xi[i], wi[i]);
    out << xi[i] << "    " << wi[i] << std::endl;
  }
  out << (k - 1) * h + h << "    " << wi[k - 1] << std::endl;
  out.close();
  std::cout << "+============================+==================+============="
               "======+========"
               "==========+\n";

  std::cout << "\nM[x] = " << mathExpectation << "\n";
  std::cout << "D[x] = " << dispersion << "\n";
  std::cout << "G[x] = " << standartDeviation << "\n\n";

  std::cout << "k = " << k << "\n";
  std::cout << "Kсв = " << freedomDegree << "\n";
  std::cout << "a = " << 0.05 << "\n\n";
  std::cout << "X2набл = " << criterionPirson << "\n";
  std::cout << "X2крит = " << 23.7 << "\n";

  if (criterionPirson < 23.7)
    std::cout
        << "Вывод: нельзя отвергать гипотезу о равномерном распределении.\n"
        << std::endl;
  else
    std::cout
        << "Вывод: можно отвергнуть гипотезу о равномерном распределении.\n"
        << std::endl;

  for (int t = 1; t <= k - 1; ++t) {
    double autocorrelation = at(k, xi, averageX, t);
    if (max < autocorrelation) {
      max = autocorrelation;
      max_i = t;
    }
    out2 << t << "    " << autocorrelation << "\n";
  }

  out2.close();

  std::cout << "max(τ) = " << max_i << std::endl;
  std::cout << "τ = " << max_i << std::endl;
  std::cout << "a(τ) = " << at(k, xi, averageX, max_i) << "\n";
  if (max_i == 1)
    std::cout << "Вывод: исследуемый ряд содержит тенденцию.\n";
  else if (max_i == k)
    std::cout << "Вывод: исследуемый ряд содержит колебания.\n";
  else
    std::cout
        << "Вывод: исследуемый ряд содержит либо тендецию, либо колебания.\n";

  std::system("gnuplot Hi2.gnu");
  std::system("gnuplot autocorrelation.gnu");
}

int main() {
  int n = 10000;

  double xMax = 0, xMin = 1;
  int k = kCoef(n);

  std::vector<double> randomNumbers;
  std::vector<double> ni(k);
  std::vector<double> xi(k);
  std::vector<double> wi(k);
  std::vector<double> ni1(k);

  int freedomDegree = k - 2 - 1;

  for (int i = 0; i < n; ++i) {
    double tempRandomNumber = generateRandomNumber();
    randomNumbers.push_back(tempRandomNumber);
  }

  std::sort(randomNumbers.begin(), randomNumbers.end());
  xMin = randomNumbers[0];
  xMax = randomNumbers[n - 1];

  // h - step
  double h = hStep(xMax, xMin, k);

  calcXiNi(xMin, xMax, n, k, h, randomNumbers, xi, ni);

  calcWi(n, k, ni, wi);

  double averageX = avrX(k, xi);
  double mathExpectation, dispersion, standartDeviation, criterionPirson;

  mathExpectation = calcMathExpectation(k, xi, wi);
  dispersion = calcDispersion(k, xi, wi, mathExpectation);
  standartDeviation = sqrt(dispersion);
  criterionPirson =
      Hi2(n, k, h, ni, ni1, xi, mathExpectation, standartDeviation);
  printInfo(n, k, h, ni, ni1, xi, wi, mathExpectation, dispersion,
            standartDeviation, criterionPirson, freedomDegree, averageX);

  return 0;
}