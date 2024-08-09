#include <iostream>
#include <z3_api.h> // IWYU pragma: keep

int main() {
  std::cout << Z3_get_full_version();
  return 0;
}
