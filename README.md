# MOFKit/ Porous material kit 开发

## 去除溶剂、金属离子

### 涉及相关的工具包
* [命令行参数解析 - cmdline](https://github.com/tanakh/cmdline)

### 编译
* C++ 11
* g++ main.cpp -o mofkit -std=c++11


### 补充
* [原子分数坐标](https://baike.baidu.com/item/%E5%8E%9F%E5%AD%90%E5%88%86%E6%95%B0%E5%9D%90%E6%A0%87/8578487?fr=aladdin)
* 原子团的划分
  * 算法实现基于`并查集(disjoint set union)`，暂未优化（路径压缩），时间复杂度分析$O(nlogh)$，空间复杂度为$O(n)$，其中`n`为节点个数，`h`为并查集构造过程中多叉树的最大高度
  * 另外`染色法(graph labeling)`也可以实现，时间复杂度$O(n^2)$，其中`n`为节点个数，空间复杂度为$O(n)$

### 问题
* 关于`cmp_solvent_known`函数的意义，与现有溶剂比较的意义？
* -s 参数的真正需求，导出溶剂？
* -m 参数的真正需求，删去溶剂？
