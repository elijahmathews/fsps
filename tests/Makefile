F90 = gfortran
F90FLAGS = -O3 -cpp -fPIC -I../src

# Filter out ALL standalone drivers (simple, autosps, lesssimple, spec_bin)
SRC_OBJS = $(filter-out ../src/simple.o ../src/autosps.o ../src/lesssimple.o ../src/spec_bin.o, $(wildcard ../src/*.o))

OBJS = test_utils.o test_driver.o

all: test_driver

test_driver: $(OBJS)
	$(F90) $(F90FLAGS) -o test_driver $(OBJS) $(SRC_OBJS)

test_utils.o: test_utils.f90
	$(F90) $(F90FLAGS) -c test_utils.f90

test_driver.o: test_driver.f90 test_utils.o
	$(F90) $(F90FLAGS) -c test_driver.f90

clean:
	rm -f *.o *.mod test_driver
