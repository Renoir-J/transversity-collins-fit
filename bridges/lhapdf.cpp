// Project-owned C ABI around the unmodified, source-locked LHAPDF library.
#include <LHAPDF/LHAPDF.h>
#include <LHAPDF/GridPDF.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>

namespace {
std::mutex setup_mutex;
std::string bound_root;
void failure(char* message, size_t length, const char* text) {
  if (message && length) std::snprintf(message, length, "%s", text);
}
bool supported(int pid) { return pid == 21 || (pid != 0 && std::abs(pid) <= 5); }
}

extern "C" const char* collins_lha_version() { return LHAPDF_VERSION; }

extern "C" void* collins_lha_open(const char* root, int policy, char* message, size_t length) {
  try {
    std::lock_guard<std::mutex> guard(setup_mutex);
    if (!root || !*root || (policy != 0 && policy != 1)) throw std::runtime_error("invalid CT10 setup");
    if (!bound_root.empty() && bound_root != root)
      throw std::runtime_error("use an isolated process for a different input directory");
    if (bound_root.empty()) {
      bound_root = root;
      LHAPDF::pathsPrepend(bound_root);
      LHAPDF::setVerbosity(0);
    }
    std::unique_ptr<LHAPDF::GridPDF> pdf(new LHAPDF::GridPDF("CT10nlo", 0));
    pdf->setInterpolator(std::string("logcubic"));
    pdf->setExtrapolator(std::string(policy == 0 ? "error" : "continuation"));
    pdf->setForcePositive(1); // Explicit CT10nlo DataVersion 4 metadata.
    return pdf.release();
  } catch (const std::exception& error) {
    failure(message, length, error.what());
    return nullptr;
  }
}

extern "C" void collins_lha_close(void* context) { delete static_cast<LHAPDF::GridPDF*>(context); }

extern "C" int collins_lha_limits(void* context, double* limits) {
  if (!context || !limits) return 1;
  auto* pdf = static_cast<LHAPDF::GridPDF*>(context);
  limits[0] = pdf->xMin(); limits[1] = pdf->xMax();
  limits[2] = pdf->q2Min(); limits[3] = pdf->q2Max();
  return 0;
}

extern "C" int collins_lha_batch(void* context, const int* pids, int nf, const double* x,
                                const double* q2, size_t count, double* output,
                                char* message, size_t length) {
  try {
    if (!context || nf <= 0) throw std::runtime_error("invalid CT10 context or flavor count");
    auto* pdf = static_cast<LHAPDF::GridPDF*>(context);
    for (size_t i = 0; i < count; ++i) {
      if (!std::isfinite(x[i]) || !std::isfinite(q2[i]) || x[i] <= 0 || x[i] > 1 || q2[i] <= 0)
        throw std::runtime_error("invalid CT10 coordinates");
      for (int f = 0; f < nf; ++f) {
        if (!supported(pids[f])) throw std::runtime_error("unsupported CT10 flavor");
        output[f + nf*i] = pdf->xfxQ2(pids[f], x[i], q2[i]);
        if (!std::isfinite(output[f + nf*i])) throw std::runtime_error("nonfinite CT10 result");
      }
    }
    return 0;
  } catch (const std::exception& error) {
    failure(message, length, error.what());
    return 1;
  }
}
