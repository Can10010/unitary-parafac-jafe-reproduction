# Unitary PARAFAC Algorithm for JAFE (Reproduction)

本仓库为 **Unitary PARAFAC 方法在联合角度–频率估计（JAFE）问题中的复现实现**，  
基于公开论文的算法思路，使用 MATLAB 进行独立实现与仿真验证。

---

## 参考论文

Lingyun Xu, Fangqing Wen, Xiaofei Zhang  
**A Novel Unitary PARAFAC Algorithm for Joint DOA and Frequency Estimation**  
IEEE Communications Letters, vol. 23, no. 4, Apr. 2019  
DOI: 10.1109/LCOMM.2019.2896593

> 本代码并非作者官方代码，而是基于论文描述的独立复现，仅用于学习与研究目的。

---

## 仓库内容

- `UPARAFAC.m`  
  基于 **前后向扩展 + 酉变换（Unitary Transform）** 的 PARAFAC 模型实现，用于联合 DOA 与载频估计。

---

## 方法简介（简要）

该实现主要包括以下步骤：

1. 构建联合空间–频率阵列信号模型  
2. 前后向扩展以增强数据结构  
3. 酉变换将复值模型转化为实值张量模型  
4. 基于 PARAFAC 的三线性分解  
5. 通过相位差与最小二乘方法估计 DOA 与频率

算法流程与论文中提出的方法保持一致，但在实现细节上可能有所差异。

---

## 参数说明（可直接修改）

在 `UPARAFAC.m` 文件开头可调整以下参数以进行不同实验：

- `M`：阵元数  
- `N`：延迟通道数  
- `tau`：延迟时间  
- `d`：阵元间距  
- `f`：信源频率  
- `theta`：DOA 角度  
- `L`：快拍数  
- `SNR`：信噪比  

代码中已给出不等式约束条件与参数选择说明，便于复现实验或自行调整。

---

## 仿真结果示例

下图展示了一组联合 DOA–频率估计结果，其中黑色星号为估计值，红色叉号为真实值：

<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/5e89d0d5-e983-4241-983a-5db7ac930997" />


（如需复现该结果，可直接运行 `UPARAFAC.m`）

---
## 依赖说明
**本仓库复现代码涉及两个外部依赖函数（均不随仓库分发），请自行获取并加入 MATLAB 路径。**
### 1) Khatri-Rao 积：`khatrirao.m`
本代码中使用了 Khatri-Rao 积运算，依赖 MATLAB Tensor Toolbox 中的函数 `khatrirao.m`。


源码及许可证请见：
http://www.tensortoolbox.org/

### 2) 复值 PARAFAC 拟合：`comfac.m`
本复现代码调用了外部函数 `comfac.m`，其版权与许可归原作者所有。

获取方式：https://ucphchemometrics.com/algorithm-for-fitting-complex-valued-parafac-cpd-model/  



---
## 声明

- 本代码仅用于**学习、研究与算法复现**
- 若用于学术研究或工程应用，请务必引用原论文
- 欢迎交流与指正实现中的问题
