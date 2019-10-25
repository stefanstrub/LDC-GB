src = orbits.cc common.cc 
obj = $(src:.cc=.o)

all: liborbits.a

liborbits.a:$(obj)
	ar rcs $@ $^

%.o: %.cc
	g++ -c $^ -o $@

clean:
	rm -f $(obj) liborbits.a
