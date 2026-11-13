/* Explicitly synthetic ABI fixture. No physics implementation or predictions. */
#include <stdint.h>
#include <stdio.h>

static int active = 0;
static int closes = 0;

int cf_hoppet_abi(void) { return 1; }
int cf_test_active(void) { return active; }
int cf_test_closes(void) { return closes; }

#define OPEN_ARGS double y, double x, double dy, int order, int layout, \
    double eps, double du, int nf, int nt, const double *tau, \
    int *n, int64_t *handle, int *status
#define IGNORE_ARGS (void)y; (void)x; (void)dy; (void)order; (void)layout; \
    (void)eps; (void)du; (void)nf; (void)nt; (void)tau

void cf_hoppet_open_batched(OPEN_ARGS) {
    IGNORE_ARGS;
    active = 1;
    *n = 2;
    *handle = 17;
    *status = 0;
    fputs("cgauss64D113.1TOO HIGH ACCURACY REQUIRED\n", stderr);
    fflush(stderr);
}

void cf_hoppet_open(OPEN_ARGS) {
    cf_hoppet_open_batched(y, x, dy, order, layout, eps, du, nf, nt, tau, n, handle, status);
}

void cf_hoppet_close(int64_t handle, int *status) {
    (void)handle;
    active = 0;
    closes++;
    *status = 0;
}

void cf_hoppet_grid(int64_t handle, int n, double *x, int *status) {
    (void)handle;
    for (int i = 0; i < n; i++) x[i] = 0.1 + 0.9 * i / (n - 1);
    *status = 0;
}

void cf_hoppet_info(int64_t handle, int n, int *nodes, int64_t *payload, int *status) {
    (void)handle;
    *nodes = n;
    *payload = 0;
    *status = 0;
}

/* The constructor must reject before any of these symbols can be called. */
void cf_hoppet_forward(void) {}
void cf_hoppet_adjoint(void) {}
void cf_hoppet_generator(void) {}
void cf_hoppet_evaluate(void) {}
void cf_hoppet_rows(void) {}
void cf_hoppet_collins(void) {}
void cf_hoppet_collins_rows(void) {}
