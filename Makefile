.POSIX:


CFLAGS = `pkg-config --cflags lattice_symmetries`

all: libplugin.a
plugin.o: plugin.c plugin.h
libplugin.a: plugin.o
	$(AR) rcs $@ $^
clean:
	rm -f plugin.o libplugin.a
