#include <fstream>
#include <iostream>
#include <random>
#include <vector>

double generateRandomNumber() {
  static std::random_device rd;
  static std::mt19937 mt(rd());
  static std::uniform_real_distribution<double> dist(0, 1);
  return dist(mt);
}

int factorial(int n) {
  if (n == 0) return 1;

  int fact = 1;

  for (int i = 1; i <= n; i++) {
    fact *= i;
  }

  return fact;
}

int CnK(int n, int k) {
  return factorial(n) / (factorial(k) * factorial(n - k));
}

std::vector<double> withoutReturn(int balls, int getBalls, int redBalls) {
  std::vector<double> p = {1, 1, 1, 1, 1, 1, 1, 1};

  for (int i = 0; i < redBalls; i++) {
    int count = i;
    double countRed = redBalls;
    double countBalls = balls;
    double exceptRed = balls - redBalls;

    for (int j = 0; j < getBalls; j++) {
      if (count-- > 0) {
        p[i] *= countRed-- / countBalls--;
      } else {
        p[i] *= exceptRed-- / countBalls--;
      }
    }

    p[i] *= CnK(getBalls, i);
  }

  return p;
}

std::vector<double> withReturn(int balls, int getBalls, int redBalls) {
  std::vector<double> p = {1, 1, 1, 1, 1, 1, 1, 1};

  for (int i = 0; i < redBalls; i++) {
    int count = i;
    double exceptRed = balls - redBalls;

    for (int j = 0; j < getBalls; j++) {
      if (count-- > 0) {
        p[i] *= (double)redBalls / balls;
      } else {
        p[i] *= exceptRed / balls;
      }
    }

    p[i] *= CnK(getBalls, i);
  }

  return p;
}

void printVec(const std::vector<double>& vec) {
  for (auto v : vec) {
    std::cout << " " << v << " ";
  }
}

void printVec(const std::vector<int>& vec) {
  for (auto v : vec) {
    std::cout << " " << v << " ";
  }
}

double checkPi(const std::vector<double>& p) {
  double sum = 0;
  for (auto i : p) {
    sum += i;
  }

  return sum;
}

std::vector<double> distributionFunction(const std::vector<double>& p,
                                         int redBalls) {
  std::vector<double> result(redBalls, 0);

  for (int i = 0; i < redBalls; i++) {
    for (int j = 0; j < i; j++) {
      result[i] += p[j];
    }
  }

  return result;
}

std::vector<int> checkpoint(const std::vector<double>& p, int redBalls, int n) {
  double random;
  std::vector<int> result(redBalls, 0);

  for (int j = 0; j < n; j++) {
    random = generateRandomNumber();
    for (int i = 0; i < redBalls; i++) {
      if (random >= p[i] && random < p[i + 1]) {
        result[i]++;
      }
    }
  }

  return result;
}

void printInFile(const std::vector<int>& vec1, const std::vector<double>& vec2,
                 const std::vector<int>& vec3, const std::vector<double>& vec4,
                 int redBalls, int n) {
  std::ofstream out("Results.txt");

  for (int i = 0; i < redBalls; i++) {
    out << i << "    " << vec1[i] << "    " << vec2[i] * n << "    " << vec3[i]
        << "    " << vec4[i] * n << std::endl;
  }

  system("gnuplot lab2.gnu");
  system("gnuplot lab2_2.gnu");
}

int main() {
  int balls = 21, getBalls = 7, redBalls = 8, n = 100000;
  std::vector<int> x = {0, 1, 2, 3, 4, 5, 6, 7};

  std::vector<double> piWOR = withoutReturn(balls, getBalls, redBalls);
  std::vector<double> piWR = withReturn(balls, getBalls, redBalls);

  std::vector<double> piWithoutReturn = distributionFunction(piWOR, redBalls);
  std::vector<double> piWithReturn = distributionFunction(piWR, redBalls);

  piWithoutReturn.push_back(1);
  piWithReturn.push_back(1);

  std::vector<int> resultWithoutReturn =
      checkpoint(piWithoutReturn, redBalls, n);
  std::vector<int> resultWithReturn = checkpoint(piWithReturn, redBalls, n);

  printInFile(resultWithoutReturn, piWOR, resultWithReturn, piWR, redBalls, n);

  std::cout << "\nPi without return {";
  printVec(piWOR);
  std::cout << "}" << std::endl;

  std::cout << "Result without return {";
  printVec(resultWithoutReturn);
  std::cout << "}\n" << std::endl;

  std::cout << "Pi with return {";
  printVec(piWR);
  std::cout << "}" << std::endl;

  std::cout << "Result with return {";
  printVec(resultWithReturn);
  std::cout << "}\n" << std::endl;

  return 0;
}