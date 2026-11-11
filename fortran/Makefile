FC = gfortran
FFLAGS = -O3 -fPIC -Wall -Wextra -Werror

all: libcollins_operator.so

libcollins_operator.so: operator_bridge.f90
	$(FC) $(FFLAGS) -shared $< -o $@

