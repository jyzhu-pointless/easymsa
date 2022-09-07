# EasyMSA

多序列比对工具包

开发者：朱瑾煜（2000012180@stu.pku.edu.cn）、叶昭辰（2000012184@stu.pku.edu.cn）

## 项目架构

```
easymsa
├── easymsa_basicdp.cpp
├── easymsa_hc.cpp
├── make.sh
├── README.md
└── static
    ├── BLOSUM62.txt
    ├── PAM250.txt
    └── file
        ├── igha1.fasta
        ├── igha2.fasta
        ├── ighd.fasta
        ├── ighg1.fasta
        ├── ighg2.fasta
        ├── ighg3.fasta
        ├── ighg4.fasta
        └── ighm.fasta
```

## 概述

我们使用两种包含动态规划的算法开发这一软件包。

软件包中共有两个 C++ 源码文件：`easymsa_basicdp.cpp`，使用直接从双序列比对推广而来的基本动态规划方法；`easymsa_hc.cpp`，使用分层聚类（hierarchical clustering）方法。

用户不需要手动编译我们的源码。只需要在 Linux/macOS 环境下，运行脚本 `make.sh`：

```shell
./make.sh
```

该脚本将自动完成编译。

`static` 文件夹中存放得分矩阵和测试序列等静态数据。

## 基本动态规划方法：`easymsa_basicdp.cpp`

### 算法原理

对于 $N$ 条序列，直接对其进行 $N$ 维动态规划，递推求其部分序列的匹配得分最大值。匹配得分为各位点的双重序列比对得分之和。

不难看出，这种方法需要的时间复杂度为 $O(\mathrm{length}^N)$，空间复杂度也为 $O(\mathrm{length}^N)$。这样的开销非常大。基于传统的动态规划算法，Carrillo 和 Lipman 提出了能限制动态规划运算空间的 Carrillo-Lipman 算法，但其复杂度并未显著降低。

文献表明，这种基本的动态规划算法只能对最多 10 条序列、最长 10000 字节长度的序列进行比对。考虑到我们的计算机性能有限，本方法只提供对最多 4 条序列、最长 500 字节长度的序列比对。

这部分代码完全为作者原创。

### 使用方法

```shell
EasyMSA -- A Light-weight Multiple Sequence Alignment Tool (Version 0.1.0)
USAGE: easymsa --type [prot/nucl] <default: prot> 
               --input [/path/to/seqA] [/path/to/seqB] [/path/to/seqC] ... 
               --matrix [pam250/blosum62/blast/trtr] <default: pam250(prot), trtr(nucl)> 
               --mode [global/local] <default: global>
               --gap-opening [integer, 0-100] <default: 5>
               --gap-extending [integer, 0-100] <default: 1>
NOTE: Input sequences should be at least 2 and no more than 4, in FASTA format, with 3 - 500 sites. 
```

其中 `--input` 为必须提供的参数，其他参数如果不提供，程序将采用默认值并显示 WARNING。错误的参数输入将终止程序并返回 ERROR。其他参数的说明如下：

- `--type`：序列的种类。prot为蛋白质，nucl为核酸。
- `--matrix`：打分矩阵殆种类。trtr为转换-颠换矩阵。
- `--mode`：比对模式。global为全局比对，local为局部比对。
- `--gap-opening` 和 `--gap-extending` 为空位罚分。

### 结果演示

对 4 条人工简单序列进行比对：

```bash
./msa --input A B C D --type nucl
```

程序返回结果如下：

```
TIP: Data loaded successfully.
WARNING: Option `matrix' is unset. Using default value `trtr'.
WARNING: Option `mode' is unset. Using default value `global'.
WARNING: Gap opening penalty is unset. Using default value 5.
WARNING: Gap extending penalty is unset. Using default value 1.
TIP: Allocating memory, please wait ...
========== RESULTS ==========
ATGCGC
A-GCGC
ATGTGC
ATGCG-

Score: 63
```

### 讨论

这种方法需要的时间复杂度、空间复杂度均为 $O(\mathrm{length}^N)$。即使采用 Carrillo-Lipman 算法（考虑到该算法步骤极其繁琐，本项目中未使用），也不会有明显降低。因此该方法并没有实用性。使用动态规划进行多序列比对被视为 NP 问题。

目前一般的多序列比对是使用类似“估算”的方法，通过牺牲一部分准确性来提高速率。传统的动态规划算法早已被淘汰。

正是考虑到这一点，本组两位成员进行了分工，分别尝试传统的和非传统的动态规划算法。

## 分层聚类方法：`easymsa_hc.cpp`

### 算法原理

多序列比对（MSA）是在双序列比对（PSA）的基础上开展的。对于待测的多条DNA/蛋白质序列，首先使用Needleman-Wunsch算法两两进行双序列比对，并将得分结果降序排列，即形成了Hierarchy。按照Hierarchy的顺序每次选取2条序列，通过与基准序列（若2条序列均不在结果集中，基准序列为结果集第一条序列；若只有1条不在结果集中，基准序列为在结果集中的另一条）进行PSA并通过添加空位的方式进行调整，使其整合进结果集里。逐渐让全部序列都进入最终的结果集，即为Cluster的过程。至此，MSA完成(1)。

### 输入数据

![image-20220620002117347](https://tva1.sinaimg.cn/large/e6c9d24egy1h3e0iidu4sj20n404it9i.jpg)

#### 序列类型（type）

​       1 – DNA；2 – 蛋白质

#### 打分矩阵（score）

用户可以选择使用内置的打分矩阵，也可以使用自定义打分矩阵。前者需要输入打分矩阵对应编号，后者需要提供自定义打分矩阵路径。

##### DNA序列比对（type = 1时）

​       1 – 等价矩阵；2 – 转换-颠换矩阵；3 – Blast矩阵

##### 蛋白质序列比对（type = 2时）

​       1 – PAM250；（适用于相似性较低序列）

​       2 – BLOSUM62（适用于相似性较高序列）

##### 用户自定义打分矩阵

输入打分矩阵路径，自定义打分矩阵要求：

①  第一行为矩阵的列名（且行名与列名顺序相同）

②  tsv格式

自定义打分矩阵可以通过R或python创建dataframe，再将dataframe以tsv格式生成文件获得。

#### 序列数量（nseq）

​       要求大于1

#### 序列文件夹路径（directory）

​       序列文件夹为包含了全部序列文件的目录。要求序列文件为FASTA格式。

​    待比对序列可以通过蛋白质/核酸序列库，比如Uniprot和NCBI搜索目标序列下载序列的FASTA文件得到。

#### 空位罚分（gap）

​       建议为正数（程序会将其转换为负分）。

#### 结果输出文件名（output）

​       可以选择将结果以文件的方式输出。 

### 输出结果

比对结果既可以直接显示在终端，也可以通过--output存入文件。终端会显示报错和提示信息（Warning和ERROR）。而结果包含比对前后的序列一一对应，空位由“-”表示。最末将显示换算后的程序运行时间（RUNNING TIME）。

![image-20220620002141717](https://tva1.sinaimg.cn/large/e6c9d24egy1h3e0iyohvxj20n409mac9.jpg)

![image-20220620002150172](https://tva1.sinaimg.cn/large/e6c9d24egy1h3e0j2a8hjj20n405cwfm.jpg)

**（输出结果展示）**

 

### 样例展示：（包含编译及运行方式）

#### 运行

​    通过给程序 `./easymsa_hc` 添加参数的方式来进行MSA。必须具备的参数为序列文件所在文件夹（--directory），其余可根据用户个人需要进行调整。

##### 获得帮助

​    命令：`./easymsa_hc`

![image-20220620002208564](https://tva1.sinaimg.cn/large/e6c9d24egy1h3e0je1pzgj20n404it9i.jpg)

##### 下载序列（以8条人类免疫球蛋白序列为例）

在Uniprot蛋白质序列网站检索igha/ighd/ighm/ighg，点击进入目标词条，点击Format选择FASTA浏览氨基酸序列，全选复制。在服务器使用mkdir file命令创建存储全部序列文件的file文件夹。进入file文件夹，用vim创建名为*.fasta（*表示用户自定义）的新文件，i键开启编辑，粘贴序列，Esc键退出编辑模式，输入”:wq”保存退出。针对每条序列重复操作得到ighA1/2.fasta、ighd.fasta、ighm.fasta、ighg1/2/3/4.fasta共计8条氨基酸序列。

**以IghM为例：**

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dza6c9usj20bd050gm0.jpg)

**（Uniprot网站检索，点击第一个词条IGHM_HUMAN）**

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dza7rcoij20bd04uglx.jpg)

**（点击Format，选择FASTA浏览序列）**

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dzaak6egj20bm01q0sk.jpg)

**（全选复制）**

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dza2l36nj20be060jrk.jpg)

**（vim创建新文件，i键编辑，粘贴，Esc键退出编辑，:wq保存退出）**

 

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dza16l51j217j01d3yu.jpg)

**（8条序列）**

 

#### 运行MSA程序

考虑进行的是蛋白质序列比对（type = 2），使用自定义的打分矩阵需要给出，访问路径（score = /www/wwwroot/easymsa/static/PAM250.txt），总共8条序列（nseq = 8），序列文件夹为/home/root2/sample/file（4.1中创建的file文件夹），空位罚分使用默认参数。

命令：./msa_hc –type 2 –score /www/wwwroot/easymsa/static/PAM250.txt –nseq 8 –directory /home/root2/sample/file

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dz9yx1ysj21270i5gu5.jpg)

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dz9xjcv8j217r0mdn71.jpg)

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dz9voms1j217o0mdqbn.jpg)

![img](https://tva1.sinaimg.cn/large/e6c9d24egy1h3dz9ukw5sj20bk04fmxd.jpg)



### 复杂度分析

令 $n$ 为序列总数量，$l$ 为序列平均长度。

##### 时间复杂度 $O(n^2l_2) + O(n^3) + O(kl^2)$

第一项为两两进行psa的时间复杂度，第二项为构建Hierarchy时排序的时间复杂度；第三项为渐进比对算法的时间复杂度。

##### 空间复杂度 $O(nl) + O(n(n-1)/2) + O(n)$

第一项为存储待比对序列的空间复杂度，第二项为存储Hierarchy的空间复杂度（容器大小为n*(n-1)/2），第三项为存储结果集的空间复杂度。

### 讨论

#### 算法运用

可以看到动态规划主要运用在 PSA 中，而在整合 PSA 的结果以进行 MSA 时采用的是渐进比对的算法。纯粹基于动态规划的 MSA 算法，如 Carrillo-Lipman 算法，虽然在当时具有开拓意义，但 $O(l^n)$ 的空间复杂度（$l$ 为序列平均长度）和接近 $O(l^n)$ 的时间复杂度（有一定程度的降低）使得动态规划现在在 MSA 中并不是最优方法，和其他算法进行结合是实现 MSA 的更有效手段。

#### 软件优化

值得进一步优化的功能为输入目标序列文件时，一方面，可以依次输入单个目标序列文件的路径，另一方面，可以直接输入一个由多个 FASTA 序列合并而成的文件。

 

**参考资料：**

1. Corpet,F. (1988) Multiple sequence alignment with hierarchical clustering. *Nucleic Acids Research*, **16**, 10881–10890.

2. Carrillo,H. and Lipman,D. (1988) The Multiple Sequence Alignment Problem in Biology. *SIAM J. Appl. Math.*, **48**, 1073–1082.
3. https://en.wikipedia.org/wiki/Multiple_sequence_alignment

4. https://blog.csdn.net/baidu_41860619/article/details/118071986