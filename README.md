PrefixSpan

==========

改进的PrefixSpan实现。

本prefixspan算法并不是原始的算法。我对它进行了改进，主要用于挖掘单条DNA中的频繁子序列。

类型：c

作用：此为将改进的prefixspan算法在mpi上并行的程序，输入为要挖掘的dna序列，输出为挖掘出的序列模式。

编译：mpicc prefixspan_mpi_v_8.c -o prefixspan_mpi_v_8

输入参数：-l 最小长度 -s 保存输出结果的文件名称 -p  要输入的序列的名称 -n 支持度 

用法：mpiexec -f host -n 10 ./prefixspan_mpi_v_8  -p ./dna/dnaY -s sava_dnaY_100_10 -l 100 -n 10

此程序要用mpicc编译，再将编译后的程序发送到各个节点，再使用mpiexec 运行，host保存着要参与运行的节点的信息，前一个-n 10表将此mpi程序在10个主机上并行，后面的-p -s -l -n 与上面的hadoop运行参数表示的意思相同
