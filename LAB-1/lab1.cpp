#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

#define IPx 1.10606015772  // intersection point
#define IPy 1

double generateRandomNumber() {
  static std::random_device rd;
  static std::mt19937 mt(rd());
  static std::uniform_real_distribution<double> dist(0, 2);
  return dist(mt);
}

inline double fx1(double x) { return sin(x); }

inline double fx2(double x) { return 2 - x; }
// inline double fx2(double x) { return M_PI / 2 + 1 - x; }

void rejectionMethod(int n) {
  double Ksi1, Ksi2, a, b;
  std::ofstream out("DistFunc.txt");

  for (int i = 0; i < n; i++) {
    a = 0;
    b = IPx;
    do {
      for (Ksi1 = generateRandomNumber(), Ksi2 = generateRandomNumber();
           (fx1(a + (b - a) * Ksi1) <= IPy * Ksi2);
           Ksi1 = generateRandomNumber(), Ksi2 = generateRandomNumber()) {
      }
    } while (0 > a + (b - a) * Ksi1 || a + (b - a) * Ksi1 > IPx);
    out << a + (b - a) * Ksi1 << "    " << Ksi2 << std::endl;
    std::cout << Ksi1 << " " << a + (b - a) * Ksi1 << "    " << Ksi2
              << std::endl;

    a = IPx;
    b = 2;
    for (Ksi1 = generateRandomNumber(), Ksi2 = generateRandomNumber();
         fx2(a + (b - a) * Ksi1) <= IPy * Ksi2;
         Ksi1 = generateRandomNumber(), Ksi2 = generateRandomNumber()) {
    }
    out << a + (b - a) * Ksi1 << "    " << Ksi2 << std::endl;
  }

  out.close();
  system("gnuplot lab1.gnu");
}

int main() {
  int n = 1;

  rejectionMethod(n);

  return 0;
}