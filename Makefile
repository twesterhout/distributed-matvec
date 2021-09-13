.POSIX:


CFLAGS = `pkg-config --cflags lattice_symmetries`

all: libplugin.a
plugin.o: plugin.c plugin.h
libplugin.a: plugin.o
	$(AR) rcs $@ $^
hello: hello.chpl states.chpl libplugin.a
	chpl -o $@ `pkg-config --cflags blas` `pkg-config --libs blas` -L. -lplugin `pkg-config --libs lattice_symmetries` plugin.h $<

clean:
	rm -f plugin.o libplugin.a
