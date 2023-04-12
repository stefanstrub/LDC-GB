src = orbits.cc common.cc 
obj = $(src:.cc=.o)
headers = $(src:.cc=.hpp) 
flags = -fPIC -I../../../common/constants/ -std=c++11

#PREFIX is environment variable, but if it is not set, then set default value
ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

all: liborbits.a

liborbits.a:$(obj)
	ar rcs $@ $^

%.o: %.cc
	g++ $(flags) -c $^ -o $@

clean:
	rm -f $(obj) liborbits.a

install: liborbits.a
	install -d $(PREFIX)/lib/
	install -m 644 $^ $(PREFIX)/lib/
	install -d $(PREFIX)/include/
	install -m 644 $(headers) $(PREFIX)/include/
