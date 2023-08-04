#include <math.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>

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

inline double hStep(double xMax, double xMin, double k) {
  return (xMax - xMin) / k;
}

inline double uiCalc(double xi, double mathExpectation,
                     double standartDeviation) {
  return (xi - mathExpectation) / standartDeviation;
}
inline double fiCalc(double u) {
  return pow(M_E, -pow(u, 2) / 2) / sqrt(2 * M_PI);
}

double avrX(int n, std::vector<double> xi) {
  double avrX = 0;
  for (int i = 0; i < n; i++) {
    avrX += xi[i];
  }

  return avrX / n;
}

double sampleVariance(int n, std::vector<double> xi, double avrX) {
  double result = 0;
  for (auto& number : xi) {
    result += pow(number - avrX, 2);
  }
  return result / n;
}

double at(int n, std::vector<double> xi, double avrX, double t) {
  double result = 0;
  double variance = sampleVariance(n, xi, avrX);
  for (int i = 0; i < n - t; i++) {
    result += ((xi[i] - avrX) * (xi[i + t] - avrX)) / (variance * (n - t));
  }

  return result;
}

int main() {
  int n = 30;
  std::vector<double> randomNumbers;
  double xMax = 0, xMin = 1;
  double k = kCoef(n) + 1;
  std::vector<double> nico(k);
  std::vector<double> xi(k);
  std::vector<double> wi(k);
  std::vector<double> niii(k);
  int freedomDegree = k - 2 - 1;
  double averageX = avrX(n, xi);
  double autocorrelation;
  // int tStep = 50;

  std::cout << "Random numbers in dist 0 and 1:"
            << "\n";
  for (int i = 0; i < n; ++i) {
    double tempRandomNumber = generateRandomNumber();
    randomNumbers.push_back(tempRandomNumber);
  }

  std::sort(randomNumbers.begin(), randomNumbers.end());
  xMin = randomNumbers[0];
  xMax = randomNumbers[n - 1];

  // h - step
  double h = hStep(xMax, xMin, k);

  double left = xMin, right = xMin + h;
  int intervalNumber = 0;

  for (auto& number : randomNumbers) {
    if (number >= right) {
      xi[intervalNumber] = (left + right) / 2;
      left = right;
      right = left + h;
      ++intervalNumber;
    }
    ++nico[intervalNumber];
  }
  xi[intervalNumber] = (left + right) / 2;

  for (int i = 0; i < k; ++i) {
    wi[i] = nico[i] / n;
  }

  std::cout << "\ncount of same intervals:\n";
  for (int i = 0; i < k; ++i) {
    std::cout << "[" << i * h << ";" << i * h + h << ") - " << nico[i] << "\n";
  }

  std::cout << "\nxi in count of same intervals:\n";
  for (int i = 0; i < k; ++i) {
    std::cout << nico[i] << " - " << xi[i] << "\n";
  }

  std::cout << "\nwi in count of same intervals:\n";
  for (int i = 0; i < k; ++i) {
    std::cout << nico[i] << " - " << wi[i] << "\n";
  }

  double mathExpectation = 0, dispersion = 0, standartDeviation = 0,
         criterionPirson = 0;
  for (int i = 0; i < k; i++) {
    mathExpectation += xi[i] * wi[i];
    dispersion += pow(xi[i], 2) * wi[i];
  }
  dispersion -= pow(mathExpectation, 2);
  standartDeviation = sqrt(dispersion);

  for (int i = 0; i < k; ++i) {
    niii[i] = ((n * h) / standartDeviation) *
              fiCalc(uiCalc(xi[i], mathExpectation, standartDeviation));
  }

  for (int i = 0; i < k; ++i) {
    criterionPirson += (nico[i] - niii[i]) / niii[i];
  }
  // evolution of variable definition:
  // iii -> someSeriousShit -> nico + nico - niii

  std::cout << "\nM[x] = " << mathExpectation << "\n";
  std::cout << "D[x] = " << dispersion << "\n";
  std::cout << "G[x] = " << standartDeviation << "\n";

  std::cout << "\nX2 = " << criterionPirson << "\n";
  std::cout << "Freedom = " << freedomDegree << "\n";

  for (int t = k - 1; t >= 1; --t) {
    autocorrelation = at(k, xi, averageX, t);
    std::cout << "\na(" << t << ") = " << autocorrelation << "\n";
  }
  std::cout << "\na(" << criterionPirson
            << ") = " << at(k, xi, averageX, criterionPirson) << "\n";

  /* std::cout <<
  "+==================+===================+===================+==="
               "===============+\n";
  std::cout
      << "|       x          |       y(x)        |       y'(x)       |        "
         "S(x)      |\n";
  std::cout
      << "+==================+===================+===================+========"
         "==========+\n";
  printf(
      "|     x = %.1lf      |  y(%.1lf)=%lf  |  y'(%.1lf)=%lf |  S(%.1lf)=%lf "
      "|\n",
      x, x, yy[0], x, yy[1], x, p1); */

  return 0;
}

/* --------------------------------------------------------------------------------------

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

inline double uiCalc(double xi, double mathExpectation,
                     double standartDeviation) {
  return (xi - mathExpectation) / standartDeviation;
}
inline double fiCalc(double u) {
  return pow(M_E, -pow(u, 2) / 2) / sqrt(2 * M_PI);
}

void calcXiNiWi(double xMin, double xMax, double n, int k, double h,
                std::vector<double>* randomNumbers, std::vector<double>* xi,
                std::vector<double>* ni, std::vector<double>* wi) {
  double left = xMin, right = xMin + h;
  int intervalNumber = 0;

  for (auto& number : *randomNumbers) {
    if (number >= right) {
      (*xi)[intervalNumber] = (left + right) / 2;
      left = right;
      right = left + h;
      ++intervalNumber;
    }
    ++(*ni)[intervalNumber];
  }
  (*xi)[intervalNumber] = (left + right) / 2;

  for (int i = 0; i < k; ++i) {
    (*wi)[i] = (*ni)[i] / n;
  }
}

double avrX(int k, std::vector<double>* xi) {
  double avrX = 0;
  for (int i = 0; i < k; i++) {
    avrX += (*xi)[i];
  }

  return avrX / k;
}

double sampleVariance(int n, std::vector<double>* xi, double avrX) {
  double result = 0;
  for (auto& number : *xi) {
    result += pow(number - avrX, 2);
  }
  return result / n;
}

double at(int n, std::vector<double>* xi, double avrX, double t) {
  double result = 0;
  double variance = sampleVariance(n, xi, avrX);
  for (int i = 0; i < n - t; i++) {
    result +=
        (((*xi)[i] - avrX) * ((*xi)[i + t] - avrX)) / (variance * (n - t));
  }

  return result;
}

double calcMathExpectation(int k, std::vector<double>* xi,
                           std::vector<double>* wi) {
  double result = 0;

  for (int i = 0; i < k; i++) {
    result += (*xi)[i] * (*wi)[i];
  }

  return result;
}

double calcDispersion(int k, std::vector<double>* xi, std::vector<double>* wi,
                      double mathExpectation) {
  double result = 0;

  for (int i = 0; i < k; i++) {
    result += pow((*xi)[i], 2) * (*wi)[i];
  }

  return result - pow(mathExpectation, 2);
}

double Hi2(int n, int k, double h, std::vector<double>* ni,
           std::vector<double>* ni1, std::vector<double>* xi,
           double mathExpectation, double standartDeviation) {
  double result = 0;

  for (int i = 0; i < k; ++i) {
    (*ni1)[i] = ((n * h) / standartDeviation) *
                fiCalc(uiCalc((*xi)[i], mathExpectation, standartDeviation));
  }

  for (int i = 0; i < k; ++i) {
    result += ((*ni)[i] - (*ni1)[i]) / (*ni1)[i];
  }

  return result;
}

void printInfo(int n, int k, double h, std::vector<double>* ni,
               std::vector<double>* ni1, std::vector<double>* xi,
               std::vector<double>* wi, double mathExpectation,
               double dispersion, double standartDeviation,
               double criterionPirson, double freedomDegree, double averageX) {
  std::ofstream out("data.txt");

  std::cout << "Random " << n << " numbers in dist 0 and 1:"
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
        i * h, i * h + h, (*ni)[i], (*xi)[i], (*wi)[i]);
    out << i * h << "    " << (*wi)[i] << std::endl;
  }
  out.close();
  std::cout << "+============================+==================+============="
               "======+========"
               "==========+\n";

  std::cout << "\nM[x] = " << mathExpectation << "\n";
  std::cout << "D[x] = " << dispersion << "\n";
  std::cout << "G[x] = " << standartDeviation << "\n";

  std::cout << "\nX2 = " << criterionPirson << "\n";
  std::cout << "Freedom = " << freedomDegree << "\n\n";

// for (int t = k - 1; t >= 1; --t) {
  //  autocorrelation = at(k, xi, averageX, t);
  //  std::cout << "\na(" << t << ") = " << autocorrelation << "\n";
 // }
std::cout << "τ = " << criterionPirson;
std::cout << "\na(τ) = " << at(k, xi, averageX, criterionPirson) << "\n";

std::system("gnuplot graph.gnu");
}

int main() {
  int n = 10000;

  double xMax = 0, xMin = 1;
  int k = kCoef(n) + 1;

  std::vector<double> randomNumbers;
  std::vector<double> ni(k);
  std::vector<double> xi(k);
  std::vector<double> wi(k);
  std::vector<double> ni1(k);

  int freedomDegree = k - 2 - 1;
  double averageX = avrX(k, &xi);
  // double autocorrelation;
  //  int tStep = 50;

  for (int i = 0; i < n; ++i) {
    double tempRandomNumber = generateRandomNumber();
    randomNumbers.push_back(tempRandomNumber);
  }

  std::sort(randomNumbers.begin(), randomNumbers.end());
  xMin = randomNumbers[0];
  xMax = randomNumbers[n - 1];

  // h - step
  double h = hStep(xMax, xMin, k);

  calcXiNiWi(xMin, xMax, n, k, h, &randomNumbers, &xi, &ni, &wi);

  double mathExpectation, dispersion, standartDeviation, criterionPirson;

  mathExpectation = calcMathExpectation(k, &xi, &wi);
  dispersion = calcDispersion(k, &xi, &wi, mathExpectation);
  standartDeviation = sqrt(dispersion);
  criterionPirson =
      Hi2(n, k, h, &ni, &ni1, &xi, mathExpectation, standartDeviation);

  printInfo(n, k, h, &ni, &ni1, &xi, &wi, mathExpectation, dispersion,
            standartDeviation, criterionPirson, freedomDegree, averageX);

  return 0;
}*/