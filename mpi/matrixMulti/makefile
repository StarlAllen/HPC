CC = mpic++
CCFLAGS = -O2 -fopenmp
LDFLAGS = -lopenblas
gemm: gemm.cpp
	${CC} ${CCFLAGS} gemm.cpp -o gemm ${LDFLAGS}
clean:
	rm gemm
