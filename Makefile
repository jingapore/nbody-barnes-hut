# wrapper to enable us to run `make` in top most directory

.PHONY: all clean nbody

all: nbody

nbody:
	@$(MAKE) -C build
	@cp build/nbody .

clean:
	@$(MAKE) -C build clean
	@rm -f nbody