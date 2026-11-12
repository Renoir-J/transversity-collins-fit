// Separate executable/reference path, using the public factory API.
#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/GridPDF.h>
#include <iomanip>
#include <iostream>
#include <memory>

int main(int argc, char** argv) {
  try {
    if (argc != 3) return 2;
    LHAPDF::pathsPrepend(argv[1]);
    LHAPDF::setVerbosity(0);
    std::unique_ptr<LHAPDF::PDF> owner(LHAPDF::mkPDF("CT10nlo", 0));
    auto* pdf = dynamic_cast<LHAPDF::GridPDF*>(owner.get());
    if (!pdf) return 3;
    pdf->setInterpolator(std::string("logcubic"));
    pdf->setExtrapolator(std::string(std::string(argv[2]) == "strict" ? "error" : "continuation"));
    pdf->setForcePositive(1);
    int count, pid;
    double x, q2;
    std::cin >> count;
    for (int i = 0; i < count; ++i) {
      if (!(std::cin >> pid >> x >> q2)) return 4;
      std::cout << "RESULT " << std::setprecision(17) << pdf->xfxQ2(pid, x, q2) << "\n";
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << error.what() << "\n";
    return 1;
  }
}
