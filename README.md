# Minimum-Weight-Spanning-Tree-Comparator
---
## 🌐 Language: [English](#english) | [中文](#中文) | [日本語](#日本語)

---

<a name="english"></a>
## English

### Description

The goal of this project is to compare the output trees of Prim's algorithm and Kruskal's algorithm on edge-weighted graphs. The program accepts an edge-weighted graph represented as an adjacency matrix and returns whether the MSTs generated by the two algorithms are identical or not. The purpose of this project is to help users understand the differences and similarities between Prim's algorithm and Kruskal's algorithm when dealing with edge-weighted graphs, and to further deepen their understanding of concepts related to the minimum spanning tree in graph theory.

#### Features

- Determines the minimum weight spanning tree.
- Compares the results of Prim's algorithm and Kruskal's algorithm.

#### Usage

**Input:** An n x n array G, of type double, representing an edge-weighted graph.

**Output:** A boolean value indicating whether the results of the two algorithms match.

**Input Example:**
'''text
7
0 1 2 0 0 0 0
1 0 3 4 0 0 0
2 3 0 0 0 0 6
0 4 0 0 4 5 0
0 0 0 4 0 0 0
0 0 0 5 0 0 0
0 0 6 0 0 0 0

---

<a name="中文"></a>
## 中文

### 项目描述

此项目的目标是比较带权边图上的普里姆算法和克鲁斯卡尔算法的输出树。该程序接受一个表示为邻接矩阵的带权边图，并返回这两种算法生成的MST是否相同。这个项目的目的是帮助用户了解Prim算法和Kruskal算法在处理边加权图时的差异和相似之处，进一步深入理解图论中与最小生成树相关的概念。

#### 功能

- 确定最小权重生成树。
- 比较普里姆算法和克鲁斯卡尔算法的结果。

#### 使用方法

**输入：** 一个 n x n 的数组G，数据类型为双精度浮点数，表示一个带权图。

**输出：** 一个布尔值，表示两种算法的结果是否匹配。

**输入示例：**
'''text
7
0 1 2 0 0 0 0
1 0 3 4 0 0 0
2 3 0 0 0 0 6
0 4 0 0 4 5 0
0 0 0 4 0 0 0
0 0 0 5 0 0 0
0 0 6 0 0 0 0


---

<a name="日本語"></a>
## 日本語

### 説明

このプロジェクトの目的は、エッジに重みが付けられたグラフ上のプリムのアルゴリズムとクラスカルのアルゴリズムの出力ツリーを比較することです。このプログラムは、隣接行列として表されるエッジの重みがついたグラフを受け入れ、2つのアルゴリズムによって生成されたMSTが同一であるかどうかを返します。このプロジェクトの目的は、ユーザーが辺の重み付きグラフを取り扱う際のPrimのアルゴリズムとKruskalのアルゴリズムの違いと類似点を理解し、グラフ理論における最小全域木に関連する概念の理解をさらに深めることです。

#### 機能

- 最小の重みの生成木を決定する。
- プリムのアルゴリズムとクラスカルのアルゴリズムの結果を比較する。

#### 使用法

**入力：** n x nの配列G、タイプはdoubleで、エッジに重みがついたグラフを表します。

**出力：** ２つのアルゴリズムの結果が一致しているかどうかを示すブール値。

**入力例：**
'''text
7
0 1 2 0 0 0 0
1 0 3 4 0 0 0
2 3 0 0 0 0 6
0 4 0 0 4 5 0
0 0 0 4 0 0 0
0 0 0 5 0 0 0
0 0 6 0 0 0 0
