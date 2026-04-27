# 1D-DVR
1D-DVR
# 1D-DVR
# 1D-DVR: One-Dimensional Discrete Variable Representation Solver

这是一个基于 Python 实现的通用一维离散变量表象（Discrete Variable Representation, DVR）求解器。该程序旨在解决一维薛定谔方程，特别适用于双原子分子或小分子振动模的本征能量和波函数计算。(写着玩的)

## 🚀 项目亮点

- **高度模块化**：采用解耦设计，底层物理引擎 (`dvr_engine.py`) 与用户运行端 (`run_morse.py`) 相互独立，可轻松适配各种一维势能函数。
- **物理基准验证**：程序已通过真实 $H_2$ 分子莫尔斯势（Morse Potential）的基准测试，零点能（ZPE）计算结果与实验值高度吻合。
- **自动格点优化**：内置 `find_grid` 算法，根据势能面极小值点和曲率自动调整物理格点范围。
- **多体系适用**：不仅适用于 Morse 势，也适用于双势阱（如 $NH_3$ 反转振动）等非简谐体系。

## 📂 目录结构

```bash
.
├── dvr_engine.py    # 核心引擎：矩阵构建、表象变换及对角化求解
├── run_morse.py     # 运行示例：定义分子参数、调用引擎并绘图
└── README.md        # 项目说明文档

```
## 安装要求
本项目依赖 Python 核心科学计算库，你可以通过以下命令快速配置环境：
```
pip install numpy scipy matplotlib
```
NumPy: 用于构建哈密顿矩阵及线性代数求解。

SciPy: 用于势能极小值搜索及特殊数学函数处理。

Matplotlib: 用于波函数与势能面的图形化展示。
