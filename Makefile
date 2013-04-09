LIBOBJS = src/pcgtr.o src/matrix.o
LIBNAME = libpcgtr.a

EXAMPLES = examples/example1 examples/example1-minimal examples/example2

CC = gcc
CFLAGS = -Wall -O3 -ffast-math -funroll-loops

all: $(LIBNAME)

$(LIBNAME): $(LIBOBJS)
	ar rcs build/$@ $(LIBOBJS)

examples: $(LIBNAME) $(EXAMPLES)

examples/example%: examples/example%.c
	$(CC) -o $@ $@.c $(CFLAGS) -Lbuild -lpcgtr -Isrc -lm

clean:
	rm -rf $(LIBOBJS) build/$(LIBNAME) $(EXAMPLES)