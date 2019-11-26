# MOFKit/ Porous material kit 开发

## 去除溶剂、金属离子

### 编译
> C++ 11
* g++ rm_mof_solvents.cpp -o rmsol -std=c++11

### 使用
* .\rmsol -i [cif_file]
* 具体使用参照参数说明

### 涉及相关的工具包
* [命令行参数解析 - cmdline](https://github.com/tanakh/cmdline)

### 补充
* [原子分数坐标](https://baike.baidu.com/item/%E5%8E%9F%E5%AD%90%E5%88%86%E6%95%B0%E5%9D%90%E6%A0%87/8578487?fr=aladdin)
* [原子分数坐标的变换矩阵](https://blog.csdn.net/aserendipity/article/details/82940016)
* 原子团的划分
  * 算法实现基于`并查集(disjoint set union)`，暂未优化（路径压缩），时间复杂度分析$O(nlogh)$，空间复杂度为$O(n)$，其中`n`为节点个数，`h`为并查集构造过程中多叉树的最大高度
  * 另外`染色法(graph labeling)`也可以实现，时间复杂度$O(n^2)$，其中`n`为节点个数，空间复杂度为$O(n)$
* 参数说明
  * `-i, --cif_in`  MOF的cif格式文件
  * `-o, --output_path`    去除结构中溶剂的结果导出路径
  * `-d, --skin_distance`   需要使用的表层距离（系数）
  * `-f, --force`   强制去除结构中的未知溶剂，默认只去除公知溶剂 
  * `-?, --help`    帮助说明

## 寻找空间群

### 编译
* g++ find_space_groups.cpp  ./include/spglib/_build/libsymspg.so -o a -I./include/spglib/src

### 设计相关的工具包
* [spglib](https://github.com/atztogo/spglib)

### 补充
* [spglib 安装与配置](https://atztogo.github.io/spglib/install.html)
* [编译引入动态链接库](https://www.cnblogs.com/dongfangchun/p/9078751.html)
* `error while loading shared libraries: libsymspg.so.1: cannot open shared object file: No such file or directory` 因为程序默认会到/lib64/libsymspg.so不在/lib64/中。因此需要增加如下命令，让程序也到指令的目录中找库
    `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./include/spglib/_build`


### 问题
* 部分api调用有问题
  * spg_niggli_reduce
  * spg_delaunay_reduce
* 

## IN CELL 功能

### 编译
* g++ in_cell.cpp -o a -std=c++11


## 问题补充
* `Window`读取文件和`Linux`读取文件的文本换行格式问题
  * 参考-[linux和windows文本换行格式问题（^M）](https://www.cnblogs.com/feer/p/9578059.html)
* 