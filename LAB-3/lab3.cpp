#include <iostream>
#include <random>
#include <string>
#include <vector>

bool widthMaxReached = false;

int generateRandomNumber(int a, int b) {
  static std::random_device rd;
  static std::mt19937 mt(rd());
  std::uniform_int_distribution<int> dist(a, b);
  return dist(mt);
}

void generateTree(int deep, int width, std::vector<std::vector<int>>& tree) {
  for (int i = 0; i < deep; i++) {
    int size = width;

    if (i != deep - 1 || widthMaxReached) {
      size = generateRandomNumber(1, width);
    }

    if (!widthMaxReached && size == width) {
      widthMaxReached = true;
    }

    tree[i] = std::vector<int>(size);

    if (i != 0) {
      for (int j = 0; j < size; j++) {
        tree[i][j] = generateRandomNumber(1, tree[i - 1].size());
      }
    } else {
      for (int j = 0; j < size; j++) {
        tree[i][j] = 0;
      }
    }
  }
}

void printTreeMap(int deep, std::vector<std::vector<int>>& tree) {
  for (int i = 0; i < deep; i++) {
    for (auto& a : tree[i]) {
      std::cout << a << " ";
    }
    std::cout << std::endl;
  }
}

int main(int argc, char** argv) {
  int deep = std::stoi(argv[1]);
  int width = std::stoi(argv[2]);
  std::vector<std::vector<int>> tree(deep);

  generateTree(deep, width, tree);
  printTreeMap(deep, tree);

  return 0;
}