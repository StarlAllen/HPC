
CXX=mpic++
#CCFLAGS=-Wall -g
# activate for compiler optimizations:
CCFLAGS=-Wall -O3
LDFLAGS=
LIBS=-pthread


linear: main.o
	$(CXX) $(LDFLAGS) -o $@ $^

%.o: %.cpp %.h
	$(CXX) $(CCFLAGS) -c $<

%.o: %.cpp
	$(CXX) $(CCFLAGS) -c $<

clean:
	rm -f *.o linear
