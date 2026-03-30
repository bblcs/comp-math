#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>
#include <ostream>

template <typename F>
void bisection(F &&f, double a, double b, double tolerance,
               size_t max_iterations) {

  if (std::signbit(f(a)) == std::signbit(b)) {
    std::cout << "sufficient condition not met" << std::endl;
    return;
  }

  double c;

  for (int i = 0; i < max_iterations; i++) {
    // std::cout << "a: " << a << " b: " << b << "\n";
    c = (a + b) / 2;
    auto fc = f(c);

    if (std::abs(a - b) < 2 * tolerance || fc == 0) {
      std::cout << "\troot: " << c << " found after " << i << " iterations"
                << std::endl;
      return;
    }

    if (std::signbit(f(a)) != std::signbit(fc)) {
      b = c;
    } else {
      a = c;
    }
  }

  std::cout << "\troot: " << c << " found after max iterations." << std::endl;
}

template <typename F, typename G, typename V>
V simple_iterations(F &&f, G &&phi, V x0, double tolerance,
                    size_t max_iterations) {
  V cur = x0;
  V next;
  for (int i = 0; i < max_iterations; i++) {
    next = phi(cur);
    if (std::abs(next - cur) < tolerance) {
      std::cout << "\troot: " << next << " found after " << i << " iterations"
                << std::endl;
      return next;
    }
    cur = next;
  }

  std::cout << "\tdid not converge, found: " << cur << std::endl;
  return cur;
}

template <typename F, typename DF, typename V>
V newton(F &&f, DF &&df, V x0, double tolerance, size_t max_iterations) {
  return simple_iterations(
      f,
      [f, df](V cur) {
        std::cout << "df: " << df(cur) << std::endl;
        return cur - (f(cur) / df(cur));
      },
      x0, tolerance, max_iterations);
}

template <typename F>
void secant(F &&f, double x0, double x1, double tolerance,
            size_t max_iterations) {
  double prev = x0;
  double cur = x1;
  double next;
  for (int i = 0; i < max_iterations; i++) {
    auto lambda = (cur - prev) / (f(cur) - f(prev));
    next = cur - lambda * f(cur);

    if (std::abs(next - cur) < tolerance) {
      std::cout << "\troot: " << next << " found after " << i << " iterations"
                << std::endl;
      return;
    }

    prev = cur;
    cur = next;
  }

  std::cout << "\tdid not converge, found: " << cur << std::endl;
}

void iter_demo(size_t iter) {
  auto tolerance = 1e-16;
  auto f = [](double x) { return std::tan(x) - x; };
  auto df = [](double x) {
    auto cos = std::cos(x);
    return (1 / (cos * cos) - 1);
  };

  std::cout << "===================================================="
            << std::endl;
  std::cout << iter << " iterations" << std::endl;
  std::cout << "===================================================="
            << std::endl;
  std::cout << std::fixed
            << std::setprecision(std::numeric_limits<double>::max_digits10);
  std::cout << "bisection" << std::endl;
  bisection(f, 4.4, 4.5, tolerance, 48);
  std::cout << "simple iterations with phi(x) = atan(x) + pi" << std::endl;
  simple_iterations(
      f, [](double x) { return std::atan(x) + M_PI; }, 4.4, tolerance, iter);

  std::cout << "newton" << std::endl;
  newton(f, df, 4.4, tolerance, iter);
  std::cout << "secant" << std::endl;
  secant(f, 4.4, 4.5, tolerance, iter);

  // TODO
  std::cout << "newton 0 -------------------------------------------"
            << std::endl;
  newton(f, df, 0.321, tolerance, iter);
}

int main(int argc, char **argv) {
  if (argv[1] == 0) {
    iter_demo(10);
    iter_demo(100);
    iter_demo(1000);
  } else {
    auto dimensions = std::stoi(argv[1]);
    auto p = [](std::complex<double> z) { return z * z * z - std::real(1); };
    auto dp = [](std::complex<double> z) { return std::real(3) * z * z; };

    std::ofstream ofs("image.ppm", std::ios::binary);
    ofs << "P6\n" << dimensions << " " << dimensions << "\n" << 255 << "\n";

    std::complex<double> roots[] = {{1.0, 0.0},
                                    {-0.5, std::sqrt(3.0) / 2.0},
                                    {-0.5, -std::sqrt(3.0) / 2.0}};
    for (int y = 0; y < dimensions; y++) {
      for (int x = 0; x < dimensions; x++) {
        std::complex<double> z(1.5 * (2.0 * x - dimensions) / dimensions,
                               1.5 * (2.0 * y - dimensions) / dimensions);
        z = newton(p, dp, z, 1e-16, 300); // TODO adapt iter number
        size_t idx = 0;
        double mindiff = std::abs(z - roots[0]);
        for (int i = 0; i < 3; i++) {
          double diff = std::abs(z - roots[i]);
          if (diff < mindiff) {
            mindiff = diff;
            idx = i;
          }
        }

        switch (idx) {
        case 0:
          ofs.put((char)255);
          ofs.put(0);
          ofs.put(0);
          break;
        case 1:
          ofs.put(0);
          ofs.put((char)255);
          ofs.put(0);
          break;
        case 2:
          ofs.put(0);
          ofs.put(0);
          ofs.put((char)255);
          break;
        }
      }
    }

    ofs.close();
    return 0;
  }
}
