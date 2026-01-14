export SPS_HOME := $(CURDIR)

all:
	$(MAKE) -C src

test: all
	$(MAKE) -C tests
	cd tests && ./test_driver

clean:
	$(MAKE) -C src clean
	$(MAKE) -C tests clean
