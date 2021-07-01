app=${1}
if [ ${app} = "gemm" ]; then
mpirun --hostfile hostfile -np ${2} ./gemm  4024 4024 4024 0
fi
if [ ${app} = "conv" ]; then
mpirun --hostfile hostfile -np ${2} ./conv 4096 4096 ${3}
fi
if [ ${app} = "pooling" ]; then
mpirun --hostfile hostfile -np ${2} ./pooling 1024 1024
fi
