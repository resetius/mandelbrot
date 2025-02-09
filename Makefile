
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Linux)
	COMPILER?=gcc
endif
ifeq ($(UNAME_S),Darwin)
	COMPILER?=gcc
endif
ifneq (,$(findstring CYGWIN,$(UNAME_S)))
	COMPILER?=x86_64-w64-mingw32-gcc
	PKG_CONFIG_LIBDIR:=/usr/x86_64-w64-mingw32/sys-root/mingw/lib/pkgconfig
	export PKG_CONFIG_LIBDIR
endif

COMPILER ?= $(CC)

All: mandelbrot_bigfloat.exe mandelbrot.exe mandelbulb.exe julia.exe

clean:
		rm -f *.o *.exe 

mandelbulb.exe: mandelbulb_gui.o mandelbulb.o
		$(COMPILER) $^ $(CFLAGS) `pkg-config --libs gtk+-3.0 epoxy` -lm -o $@

%.exe: %.o Makefile
		$(COMPILER) $< $(CFLAGS) `pkg-config --libs gtk+-3.0 epoxy` -lm -o $@

%.o: %.c Makefile
		$(COMPILER) -g -Wall $(CFLAGS) `pkg-config --cflags gtk+-3.0 epoxy` -c $< -o $@

%.o: %.cpp Makefile
		$(COMPILER) -g -Wall -std=c++20 $(CFLAGS) `pkg-config --cflags gtk+-3.0 epoxy` -c $< -o $@

