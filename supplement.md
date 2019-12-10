# MOFKit/ Porous material kit 开发

---
---

## 去除溶剂、金属离子

### 编译
> C++ 11
* g++ rm_mof_solvents.cpp -o a -std=c++11

### 相关工具包
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

--- 

## 寻找空间群

### 编译
* g++ find_space_groups.cpp  ./include/spglib/_build/libsymspg.so -o a -I./include/spglib/src

### 相关工具包
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

---

## IN-CELL

### 编译
* g++ in_cell.cpp -o a -std=c++11

---

## ICSD 分类与去重


### 编译

### 思路
* 逐一读入输入文件cif文件，按照：component/element type/space group/（元素总数/元素类型/空间群编号）建立相应的文件夹进行归类
* 如果发现space group文件夹下有cif文件，则进行去重，按照相似性法则，如果相似，根据文献发表时间顺序，保存最新发表的结构；如果不相似，则全部保存

### 相似性计算
**仿射对应（affine mapping）**
基本原理是将两个晶体结构中的原子点阵互相对应起来，首先把需要对比的两个晶体结构调节为同一原子数密度，然后找到一个仿射转换 T 能将 A 结构中的每个原子都和 B 结构中的每个原子唯一对应上。如果两个结构可以建立起这样的仿射对应关系，则认为这两个结构是相似的并且属于同一晶体结构原型。

**配位特征函数（Coordination Characterization Function，CCF）**
一个稳定且高效的结构表征方法 CCF，通过计算原子间距并使用归一化的高斯函数将其展宽并累加得到光滑连续的 CCF。使用关联系数公式可以计算出不同结构的 CCF 的相似性，最终可以得到一个定量的数值标定两结构间的差异程度。



---

## CSD 分类

> 剔除含金属有机分子和`disorder`分子

### 编译
* g++ CSD_classify.cpp -o a -std=c++11


---

## 格式转换

> 将 cif 文件转换成特定的格式

### 编译
* g++ format_conversion.cpp -o a -std=c++11

### 说明
* 参数说明
  * `-i, --cif_in`  MOF的cif格式文件
  * `-o, --output_path`    去除结构中溶剂的结果导出路径
  * `-d, --skin_distance`   需要使用的表层距离（系数）
  * `-f, --force`   强制去除结构中的未知溶剂，默认只去除公知溶剂 
  * `-?, --help`    帮助说明


---

## 问题补充
* `Window`读取文件和`Linux`读取文件的文本换行格式问题
  * 参考-[linux和windows文本换行格式问题（^M）](https://www.cnblogs.com/feer/p/9578059.html)

