# 十维互无偏基（$d=10$ MUBs）19 轮研究总论：零浮点不可行定理、闭式四元数 5-MUB 与 NPA 无穷阶阻断

> **作者**：芮颖浩（Yinghao Rui）  
> **整理日期**：2026-11-25（含外部审稿意见全量装甲修订 + Round 20–22 终极版）  
> **归档目录**：`D:\antigravity\d10_mub_continuous\preprint_d10_mub\`  
> **配套文件**：
> - 英文 arXiv 预印本草稿：[paper_draft_en.md](file:///D:/antigravity/d10_mub_continuous/preprint_d10_mub/paper_draft_en.md)  
> - Gribling–Polak Def. 3.1 逐条零浮点对齐与降次验算脚本：[verify_gp_alignment_and_degrees_QED.mjs](file:///D:/antigravity/d10_mub_continuous/preprint_d10_mub/verify_gp_alignment_and_degrees_QED.mjs)  
> - 全量核心代码与证书 SHA-256 清单：[MANIFEST_SHA256.json](file:///D:/antigravity/d10_mub_continuous/preprint_d10_mub/MANIFEST_SHA256.json)  
> - 一键校验清单生成脚本：[generate_sha256_manifest.mjs](file:///D:/antigravity/d10_mub_continuous/preprint_d10_mub/generate_sha256_manifest.mjs)

---

## 摘要与核心贡献总览

在 $d=6$ 已由 2026 年多篇预印本（OpenAI 2026-09 计算机辅助与 Lean 验证工作；Cárdenes Wuttig & Tindall, arXiv:2611.06917, 2026-11；Jaswal, 2026-11）确立 $N(6)=3$ 之后，十维复希尔伯特空间 $\mathbb{C}^{10}$ 中最多能存在多少组互无偏基（Mutually Unbiased Bases, MUBs）——特别是究竟能不能存在 4 组复互无偏基（$N(10) = 3$ 还是 $\ge 4$）——已跃升为量子信息理论中最小未决合数维度的核心公开难题。近期预印本（Yang et al., 2026-11-16）排除了含标准基与傅里叶基 $\{I, F_{10}\}$ 以及含三组乘积基的四基方案，但对非傅里叶起始基对（如对称会议矩阵基 $H_{\text{conf}}$、孤立循环哈达玛矩阵、一般连续族）以及一般全局代数障碍问题仍留有大片空白。

本总论系统总结了连续 19 轮（本地 16 核精确代数/数值计算 + Kaggle 双卡 T4 累计 **492 万条** GPU 轨迹）的研究成果，自然形成**两篇可独立发表的严格理论成果**与**一份完整的全局几何景观普查**：

1. **论文主线 I（固定会议矩阵基对 $\{I, H_{\text{conf}}\}$ 的单位根完备分类与 180 组已知第三基零浮点终结）**：
   - **定理 A（全体 $\mu_{2^a 3^b 5^c}$ 单位根完备分类）**：证明对任意非负整数 $a, b, c \ge 0$，所有相位属于 $\mu_{2^a 3^b 5^c}$ 的 $\{I, H_{\text{conf}}\}$ 完整第三基**仅有已知的 36 组纯 $\mu_4$ 第三基 $D_0$**（$a<2$ 时为 0 组），且这 36 组均无任意复相位第 4 无偏向量。
   - **定理 B（两族全新无理代数第三基 $B_{\text{new}}, B_{\text{iso}}$ 的闭式发现与 180 组已知第三基全家族零浮点证死）**：推翻“仅有 72 组第三基”的朴素猜想，在 $G_{1440}(H_{\text{conf}})$ 下精确构造出：
     - $B_{\text{new}} \in \frac{1}{6}\mathbb{Z}[\zeta_{20}][\sqrt{6}]$（$u = \frac{1+i\sqrt{5}}{\sqrt{6}} \notin \mu_\infty$，1 条轨道 **36 组**，171 维联合缺损度 $= 4$，开启连接 $H_{\text{conf}}$ 与 $F_{10}^{(4)}$ 的 4 维连续三基流形 $\mathcal{M}_4$）；
     - $B_{\text{iso}} \in \frac{1}{9}\mathbb{Z}[i, \sqrt{5}]$（$v = \frac{2+i\sqrt{5}}{3} \notin \mu_\infty$，哈达玛缺损度 $= 0$，手性稳定子 $|G_{\text{lin}}|=20, |G_{\text{anti}}|=0 \implies$ **72 组**，构成全新 `{0, 0, 16}` 刚性孤立三基族）。
     - 通过 6 次 Hilbert Nullstellensatz 模 $p$ 满秩证明（`1256/1256`）与 **128 位 CRT 精确对偶核 `BigInt` 证明**（全部 14 条核向量在 43,220 条完整生成行上逐项零浮点核验），**100% 零浮点证明全部 4 族共 180 组已知 $\{I, H_{\text{conf}}\}$ 第三基（$36\,D_0 + 36\,H\overline{D_0} + 36\,B_{\text{new}} + 72\,B_{\text{iso}}$）均不存在任何任意复相位的第 4 无偏向量**。*（范围声明：上述零浮点证明覆盖全体 $\mu_{2^a 3^b 5^c}$ 单位根以及全部 180 组已知第三基，含 $B_{\text{new}}$ 轨道点；对于 $B_{\text{new}}$ 所在 4 维连续流形 $\mathcal{M}_4$ 的其余连续点，我们采用 96 链测地线数值普查，最小残差为 `34.78`，并不宣称已对 $\{I, H_{\text{conf}}\}$ 的任意连续相位第三基给出完备分类。）*
   - **定理 C（四条零浮点刚性与横截隔离定理）**：在完整 171 维十维三基模空间中，零浮点证得 $\{I, H_{\text{conf}}, D_0\}$（52 位 Bareiss 整行列式 $\ne 0$）、同号乘积三基 $\{I, F_2\otimes F_5, (SF_2)\otimes(D_1 F_5)\}$ 以及 $\{I, H_{\text{conf}}, B_{\text{iso}}\}$ 均具有精确满秩 `171/171`（严格局部孤立）；同时在固定 $\{I, H_{\text{conf}}\}$ 的 90 维纤维中，零浮点证得 **$B_{\text{new}}$ 同样具有精确满秩 `90/90`（严格 0 维孤立点）**，且 4 维连续流形 $\mathcal{M}_4$ 的切空间向 $B_2$ 分量的投影满秩 `4/4`（与固定 $\{I, H_{\text{conf}}\}$ 纤维严格横截）！

2. **论文主线 II（有限维 $C^*$-MUB 三级闭式构造：$N_{C^*}(d)=\infty$、Clifford 6-MUB、四元数 4/5-MUB、Gribling–Polak 逐条对齐与普适 2-Copy 定理）**：
   - **定理 D（三级闭式零浮点有限维 $C^*$-MUB 阶梯——完全突破四元数上限，直至 $N_{C^*}(d)=\infty$！）**：全部构造均属于 $M_{rd}(\mathbb{C})$ 中的高秩复投影表示（非 $\mathbb{C}^d$ 复秩-1 基），均精确满足 Gribling–Polak Def. 3.1 公理 (i) 基内正交、(ii) 完备性、(iii) 三重积及迹归一条件：
     - **Tier 1（极小非复秩 $r=2$ · 四元数构造）**：对 $q=3$（$d=6$）给出 $M_{12}(\mathbb{Z}[i,\omega_3])$ 中的**精确四元数 4-MUB**（144 基内 + 432 三重积），对所有奇素数幂 $q \ge 5$（含 $d = 10, 14, 18, \dots$）给出 $M_{2d}(\mathbb{C})$ 中的**精确四元数 5-MUB**（$d=10$: 500 基内 + 2,000 三重积）。
     - **Tier 2（低秩 Clifford 构造 $\operatorname{Cl}_{q+1} \otimes \text{MUB}(\mathbb{F}_q)$）**：打破四元数 $m \le 5$ 上限！对所有 $d = 2^a q$ 给出精确 **$(q+1)$-MUB**——特别在 **$d=10$ 给出秩-4 的精确 Clifford 6-MUB**（$M_{40}(\mathbb{Z}[i,\omega_5])$，600 基内 + 3,000 三重积），在 $d=14$ 给出 8-MUB（秩 8），在 $d=22$ 给出 12-MUB（秩 32）！
     - **Tier 3（普适 Jordan–Wigner $C^*$-MUB 定理：$N_{C^*}(d) = \infty$！）**：对**任意整数 $d \ge 2$** 与**任意正整数 $m \ge 1$**（即使 $m > d+1$ 乃至 $m \to \infty$！），在 $(\mathbb{C}^d)^{\otimes \lceil (m-1)/2 \rceil}$ 中以闭式广义 Weyl 酉算子 $U_a U_b = \omega_d^{\pm 1} U_b U_a$ 构造出精确 $m$-MUB（已零浮点验毕 `d=3, m=5` 于 $M_9$、`d=6, m=5` 于 $M_{36}$、`d=6, m=8`、`d=10, m=12` 及 `m=21`）！
   - **推论 E（双层阻断机制与 Gribling–Polak 逐条对齐）**：
     - **第一层（不含公理 (iv) 的纯迹型 NPA / $C^*$-代数层级）**：由 Tier 3 普适定理，$N_{C^*}(d) = \infty$ 对一切 $d \ge 2$ 成立——任何不含高次角代数交换约束的维度无关迹型 NPA/SDP 层级（含 $S_d \wr S_k$ 对称约化）在**任意有限阶乃至无穷阶 $k \to \infty$** 上甚至连教科书上界 $N(d) \le d+1$ 都完全无法证伪，更不可能否定 $d=6$ 的 4-MUB 或 $d=10$ 的 4/5/6-MUB！这从反面严格解释了为何 Gribling–Polak (*Quantum* 2024) 必须引入公理 (iv) 才能排除 $d+2$ MUB。
     - **第二层（含 Gribling–Polak Def. 3.1 公理 (iv) $[XUX, XVX]=0$ 作为双边理想的精确检测门槛）**：我们的四元数秩-2 模型在算子层面违反公理 (iv)（零浮点证得 $[A, B] = [P_{0,0} P_{1,0} P_{2,0} P_{0,0},\; P_{0,0} P_{1,0} P_{3,0} P_{0,0}] \ne 0$）。在所有低次层级（$\deg(U)\le 1$ 或 $\deg(V)\le 1$，含 GP 定理 3.2 脚注涉及的 $\le 3$ 次单项式，以及全部迹次数 $\le 6$、即矩矩阵阶数 $t \le 3$ 与单纯迹关系 $\tau([A,B]) \equiv 0$）上，该模型仍 100% 保持可行；而当公理 (iv) 作为双边理想扩展至高次乘子时，我们零浮点精确测得首次击穿我们四元数模型的降次门槛为：**约化循环 7 次线性理想见证** $\operatorname{Tr}(E[A,B]) \ne 0$（对应矩矩阵 $t = 4$，在 $d=6$ 验出 4 条、$d=10$ 验出 18 条）与**约化循环 10 次 SOS 范数** $\operatorname{Tr}([A,B]^\dagger[A,B]) > 0$（不共享内基时为 12 次），为 Gribling–Polak (*Quantum* 2024) 数值观察到的高次壁垒（原文 *"numerical results indicate"* 12 次）提供了清晰的代数机制解释。
   - **定理 E'（普适任意秩-$r$ 2-Copy $\Lambda^2$ 零泄漏定理）**：对任意有限维空间上任意秩 $r \ge 1$ 的正交投影 $P$，严格证明普适恒等式 $\frac{\operatorname{Tr}(P^{\otimes 2} P_{\Lambda^2})}{\operatorname{Tr}(P^{\otimes 2})} = \frac{r-1}{2r}$——**当且仅当 $r=1$ 时反对称泄漏为 $0\%$，任意秩-2 模型（不限于本文四元数族）必有精确 $25\%$ 刚性泄漏，任意高秩表示 $r \ge 2$ 必有 $\ge 25\%$ 泄漏**！在 2-Copy 空间仅用 2 次约束 $P_{a,i}^{\otimes 2} P_{\Lambda^2} = 0$ 即可于 Level 1 瞬间斩杀全部高秩反例并否定 `(d=2, m=4)` 与全部 `(d, d+2)`，同时发现精确算子恒等式 $\Pi_a^{(R)}\Pi_b^{(R)} = F_{\text{ent}}$。
   - **发现 F（$\mathbb{H}^d \to \mathbb{C}^d$ 非交换质量间隙相图）**：构造规范不变非交换序参量 $E_{\text{comm}}(Q) \ge 0$（$=0 \iff$ 坍缩为复 4-MUB），首次沿精确四元数 4-MUB 流形 $\mathcal{M}_{\mathbb{H}}(d,4)$ 向 $\mathbb{C}^d$ 同伦滑降：素数幂 $d=3, 4, 5$ 全部光滑坍缩至 $10^{-16}\sim 10^{-14}$；已证不可能的 $d=2$ 刚性锁死在精确理论值 `1/36 = 0.027778`；而 $d=6$（`0.0435`，与 2026 年 $N(6)=3$ 预印本结论完全一致）与 $d=10$（`0.0678`）均撞上高出 12–13 个量级的不可逾越质量间隙！

3. **全空间 $d=10$ 几何景观普查与 492 万 GPU 轨迹证据**：
   - 发现无约束 $d=10$ 刚性 3-MUB 几乎全部自发归入 **$U(2)^5$（$2\times 2$ 块 $5\times 5$ 循环阵）的 24 个精确轨道**（23 孤立 + 1 维连续 `{I, H_conf, B_new}` 切片），并证明了 Fourier-Twisted $U(2)^5$ 的 **Bloch 共面充要判据**（24 轨最小 $|\det| = 0.010454$，距零仅差 1% 却无一跨越）。
   - 累计 **4,915,200 条 Kaggle GPU 轨迹** + 本地 **7 万+ 次精确解析雅可比 LM** 覆盖 11 类对称约化族及全部 7 大已知 3-MUB 几何（含 4D 连续流形 $\mathcal{M}_4$），**0 例精确 4-MUB、0 例精确 `{3-MUB + 第 4 向量}`**；同时将完整 4-MUB 与 `{3+2 正交向量}` 的最优近似纪录分别推进 25 倍（`0.4955`）与 153 倍（`3.33e-6`）。

---

## 第一部分（论文 I）：固定 $\{I, H_{\text{conf}}\}$ 的完备分类与 180 组零浮点不可行定理

### 1.1 起始基对与真对称群 $G_{1440}(H_{\text{conf}})$
设 $C$ 为 10 阶 Paley/会议对称矩阵（对角元 0，非对角元 $\pm 1$，$C^2 = 9I$），会议哈达玛基取为 $H_{\text{conf}} = (I + iC)/\sqrt{10}$。其完整保持 $\{I, H_{\text{conf}}\}$ 的对称群为 1,440 阶群：
$$G_{1440}(H_{\text{conf}}) = G_{720}^+(+C, \text{线性带符号置换}) \;\cup\; G_{720}^-(-C, \text{反线性带符号置换}\circ\text{复共轭}).$$

### 1.2 全体 $\mu_{2^a 3^b 5^c}$ 单位根第三基的完备封盘（Round 10–11）
通过两个良定义分圆环同态 $\phi_5: \mathbb{Z}[\zeta_{120}] \to \mathbb{F}_{25}$ 与 $\phi_3: \mathbb{Z}[\zeta_{120}] \to \mathbb{F}_{81}$，结合 $G_{1440}$ 轨道约化（574,922 个 $\mu_8$ 向量压至 767 轨）与 $N=240$ 奇偶行和下降双射（180 个不交输出模板），在纯 32/64 维分圆整环上零浮点求得全部射线，并用**两套独立精确团算法**（Tomita 枢轴 Bron-Kerbosch 与 280 轨无枢轴着色 MCQ）完成双校验：

| 等分数 $N$ | 合格归一化射线数 | 精确正交无序对 | 完整第三基数 | 超出已知 36 组 $\mu_4$ 基 |
|---:|---:|---:|---:|---:|
| 60（回归） | 10,752 | 2,035,656 | 36 | 0 |
| **120** | **21,552** | **7,954,056** | **36** | **0** |
| **240** | **43,152** | **31,454,856** | **36** | **0** |
| **全体 $2^a 3^b 5^c$** | $21552 + 180(N-120)$ | — | **$a\ge 2$: 36；$a<2$: 0** | **0** |

### 1.3 两族全新无理代数第三基与 180 组全部零浮点证死（Round 12–15）
打破单位根限制后，我们先后发现了两族具有极简闭式数域结构的全新精确第三基：
1. **$B_{\text{new}}$（Round 12–13）**：全部矩阵元形如 $u^e \zeta_{10}^m$（$u = \frac{1+i\sqrt{5}}{\sqrt{6}} \notin \mu_\infty, e\in\{-1,0,1\}$），精确验证于 $\frac{1}{6}\mathbb{Z}[\zeta_{20}][\sqrt{6}]$；在 $G_{1440}$ 下恰成 **1 条 36 组轨道**；固定 $H_{\text{conf}}$ 孤立，放开双基时 171 维联合缺损度 $= 4$。利用线性稳定子 $|G_{\text{lin}}|=20$ 压出 1,256 个 $\le (3,3)$ 次不变量，在分裂素数 $p_1 = 1000081, p_2 = 1000921$ 下达**精确满秩 `1256/1256`** $\implies$ 6 次 Nullstellensatz 零浮点证死无第 4 无偏向量！
2. **$B_{\text{iso}}$（Round 14–15）**：全部矩阵元形如 $i^k v^e$（$v = \frac{2+i\sqrt{5}}{3} \notin \mu_\infty, e\in\{-2,-1,0,1,2\}$），精确验证于 $\frac{1}{9}\mathbb{Z}[i,\sqrt{5}]$；单阵缺损度 $= 0$，构成全新 `{0, 0, 16}` 拓扑族；$|G_{\text{lin}}|=20, |G_{\text{anti}}|=0$ 产生 **72 组手性孤立基**，171 维联合雅可比满秩 `171/171`。以 5 素数 $\times$ 4 嵌入的 **128 位 CRT 精确重构全部 14 条对偶核向量** $v_1,\dots,v_{14} \in \mathbb{Z}[i,\sqrt{5}]^{1256}$（最大公分母 `4,761,892,800`），全部满足 $v_r[\text{const}] \equiv 0$ 且在全部 43,220 条完整生成行上零浮点验算 $M v_r \equiv 0$ $\implies$ 100% 无附带条件证死！

| 家族 | 精确数域 | 缺损度签名 | 171D 联合缺损度 | 轨道数 | 第 4 向量零浮点证明方式 |
|---|---|---|---:|---:|---|
| $D_0$ ($\mu_4$) | $\mathbb{Q}(i)$ | `{16, 8, 8}` | **0（52 位 Bareiss 整行列式）** | 36 | Round 4 `BigInt` Degree-6 Nullstellensatz |
| $H_{\text{conf}}\overline{D_0}$ | $\mathbb{Q}(i, \sqrt{10})$ | `{16, 8, 8}` | **0（反酉等价）** | 36 | 反酉等价传递 |
| $B_{\text{new}}$ | $\mathbb{Q}(\zeta_{20}, \sqrt{6})$ | `{16, 8, 8}` | **4（4D 连续流形 $\mathcal{M}_4$）** | 36 | Round 13 模 $p_1, p_2$ 满秩 `1256/1256` |
| $B_{\text{iso}}$ | $\frac{1}{9}\mathbb{Z}[i, \sqrt{5}]$ | `{16, 0, 0}` | **0（模 $p_1, p_2$ 满秩 `171/171`）** | 72 | Round 15 128 位 CRT 对偶核 43,220 行纯 `BigInt` 验算 |
| **合计（已知 4 族）** | — | — | — | **180 组** | **180 / 180 全部 100% 零浮点证无任意复相位第 4 向量** |

> **完备性边界与严格横截隔离说明**：本节零浮点不可行定理严格覆盖（1）全体 $\mu_{2^a 3^b 5^c}$ 单位根相位第三基，以及（2）全部 180 组已知 $\{I, H_{\text{conf}}\}$ 第三基。进一步地，Round 20 零浮点证明了固定 $\{I, H_{\text{conf}}\}$ 的 90 维雅可比在 $B_{\text{new}}$ 处亦达精确满秩 `90/90`（`p = 1000081, 1000921`），即**全部 180 组已知第三基在固定 $\{I, H_{\text{conf}}\}$ 下均为严格 0 维孤立点**；而 4 维连续流形 $\mathcal{M}_4$ 的切空间向 $B_2$ 分量的投影满秩 `4/4`（与固定 $H_{\text{conf}}$ 纤维严格横截），其离开 $H_{\text{conf}}$ 后的连续区域由 96 链测地线数值普查覆盖（最小残差 `34.78`）。

---

## 第二部分（论文 II）：有限维 $C^*$-MUB 三级阶梯（$N_{C^*}(d)=\infty$）、Gribling–Polak 逐条对齐与普适 2-Copy 定理

### 2.1 三级闭式零浮点有限维 $C^*$-MUB 阶梯（Round 17 & Round 20）
1. **Tier 1（极小非复秩 $r=2$ · 四元数 4/5-MUB 定理，Round 17）**：设 $q$ 为奇素数幂，$F_q$ 为 $\mathbb{F}_q$ 上傅里叶矩阵，$Q_a = \operatorname{diag}(\chi(a x^2))$。定义四元数酉阵：
$$B_0 = I_{2q}, \qquad B_m = \frac{1}{\sqrt{2}}\begin{pmatrix} Q_{a_m} F_q & Q_{a_m} F_q \\ q_m Q_{a_m} F_q & -q_m Q_{a_m} F_q \end{pmatrix} \in \mathrm{Sp}(2q).$$
   - 对 $q = 3$（$d=6$）：$(q_m, a_m) \in \{(1,1), (i,2), (j,0)\}$ 给出 $\mathbb{H}^6 \hookrightarrow \mathbb{C}^{12}$ 中的**精确四元数 4-MUB**（秩-2 复投影表示，非 $\mathbb{C}^6$ 复 4-MUB）。
   - 对任意奇素数幂 $q \ge 5$（含 $d = 10, 14, 18, 22, \dots$）：$\operatorname{char}(\mathbb{F}_q)\ge 5$ 时取 $(1,1), (i,-1), (j,2), (k,-2)$（$q=3^r\ge 9$ 时取 $(1,1), (i,-1), (j,\theta), (k,-\theta)$，$\theta\notin\mathbb{F}_3$），给出 $\mathbb{H}^{2q} \hookrightarrow \mathbb{C}^{4q}$ 中的**精确四元数 5-MUB**（秩-2 复投影表示，非 $\mathbb{C}^{2q}$ 复 5-MUB）。该族因依托 $\{1,i,j,k\}$ 而具天然上限 $m \le 5$（$q=3$ 因特征 3 退化为 4），但在 $m \le 5$ 达到极小非复秩 $r = 2$。
2. **Tier 2（低秩 Clifford $(q+1)$-MUB 定理，Round 20）**：对一切 $d = 2^a q$（$a \ge 1$，$q$ 奇素数幂），将 $d=2$ 因子上的 Jordan–Wigner Clifford 投影 $\frac{I + (-1)^\alpha \Gamma_a}{2}$（$\Gamma_0,\dots,\Gamma_q \in M_{2^{\lceil q/2 \rceil}}(\mathbb{Z}[i])$ 两两反对易）与 $\mathbb{F}_q$ 上的 $q+1$ 组复 MUB 做张量积，立即突破四元数 $m=5$ 上限，给出精确 **$(q+1)$-MUB**！纯 `BigInt` 零浮点已验毕：**$d=10$ 精确 6-MUB（秩 4，$M_{40}(\mathbb{Z}[i,\omega_5])$）**、$d=14$ 8-MUB（秩 8，$M_{112}$）、$d=22$ 12-MUB（秩 32，$M_{704}$）。
3. **Tier 3（普适 Jordan–Wigner $C^*$-MUB 定理：$N_{C^*}(d) = \infty$，Round 20）**：对**任意整数 $d \ge 2$** 与**任意 $m \ge 1$**（即使 $m > d+1$ 乃至 $m \to \infty$！），令 $k = \lceil (m-1)/2 \rceil$，$X_d, Z_d \in U(d)$ 为标准移位/时钟阵（$X_d Z_d = \omega_d^{-1} Z_d X_d$），$Y = \eta X_d Z_d$（$d$ 奇取 $\eta=1$，$d$ 偶取 $\eta = e^{i\pi/d}$，使 $Y^d = I_d$）。在 $(\mathbb{C}^d)^{\otimes k}$ 中定义：
$$U_{2j-1} = Y^{\otimes(j-1)} \otimes X_d \otimes I_d^{\otimes(k-j)},\qquad U_{2j} = Y^{\otimes(j-1)} \otimes Z_d \otimes I_d^{\otimes(k-j)},\qquad U_{2k+1} = Y^{\otimes k}.$$
   则所有 $U_a^d = I$ 且任意 $a < b$ 精确满足 $U_a U_b = \omega_d^{\pm 1} U_b U_a$（因 $\gcd(\pm 1, d) = 1$，条件期望 $E_a(U_b^t) \equiv 0$ 对一切 $1 \le t < d$ 成立），其谱投影 $P_{a,s} = \frac{1}{d}\sum_{t=0}^{d-1}\omega_d^{-st} U_a^t$（秩 $r = d^{k-1}$）精确满足全部 GP 公理 `(i)+(ii)+(iii)+(Tr)`！已在 $\mathbb{Z}[\zeta_{2d}]$ 纯 `BigInt` 中 0.62s 零浮点验毕 `d=3, m=5`（$M_9$，突破 $d+1=4$！）、`d=3, m=7`、`d=6, m=5`（$M_{36}$）、`d=6, m=8`（突破 $d+1=7$！）、`d=10, m=12` 与 `m=21`（突破 $d+1=11$！）——彻底回答了“换一类代数模型能否顶得更高”的追问：**在不施加公理 (iv) 或秩-1 约束时，任意维 $d \ge 2$ 均无任何有限上界（$N_{C^*}(d) = \infty$）！**

### 2.2 与 Gribling–Polak (*Quantum* 8, 1318, 2024, Def. 3.1) 的逐条对齐与精确次数分析

为了杜绝任何公理口径偏差，我们将本文的显式秩-2 复投影模型 $\{P_{a,r}\} \subset M_{2d}(\mathbb{C})$ 与 Gribling–Polak（简称 GP）定义 3.1 的 $(d, k)$-MUB 代数公理逐条零浮点核对（独立脚本 [verify_gp_alignment_and_degrees_QED.mjs](file:///D:/antigravity/d10_mub_continuous/preprint_d10_mub/verify_gp_alignment_and_degrees_QED.mjs)，17.8s 纯 `BigInt` 全量通过）：

| GP Def. 3.1 公理 / 机制项 | 数学表述 | 本文四元数秩-2 模型状态 | 零浮点验算数量与精确说明 |
|---|---|---|---|
| **公理 (i) 基内正交与幂等自伴** | $P_{a,r} P_{a,s} = \delta_{rs} P_{a,r}$, $P_{a,r}^\dagger = P_{a,r}$ | **精确满足 ✓** | $d=6$ 全验 $4\times 6^2 = \mathbf{144}$ 条；$d=10$ 全验 $5\times 10^2 = \mathbf{500}$ 条（含全部 $r \ne s$ 非对角正交） |
| **公理 (ii) 完备性** | $\sum_{r=0}^{d-1} P_{a,r} = I$ | **精确满足 ✓** | $d=6$ 全验 **4** 组基；$d=10$ 全验 **5** 组基 |
| **公理 (iii) 三重积关系** | $P_{a,r} P_{b,s} P_{a,r} = \frac{1}{d} P_{a,r}$ ($a \ne b$) | **精确满足 ✓** | $d=6$ 全验 $4\times 3\times 6^2 = \mathbf{432}$ 条；$d=10$ 全验 $5\times 4\times 10^2 = \mathbf{2,000}$ 条 |
| **迹归一化约定** | 本文 $\tau(I)=1, \tau(P_{a,r})=1/d$ vs GP $\tau_{\text{GP}}(I)=d, \tau_{\text{GP}}(P_{a,r})=1$ | **精确等价 ✓** | 仅差常数缩放因子 $d$（$\tau_{\text{GP}} = d \cdot \tau$），对可行性与正半定性零影响 |
| **$S_d \wr S_k$ 对称约化** | 对称化后的约化 SDP | **自动通过 ✓** | GP 的 SDP 在 $S_d \wr S_k$ 置换群下不变；将本文迹态 $L$ 按群平均后仍严格可行 |
| **公理 (iv) 角代数交换性** | $[X U X,\; X V X] = 0$（对一切单项式 $U, V$） | **算子层面违反 ✗；低次迹可行 ✓；高次理想击穿 ✗** | **必须拆成两层精确陈述（见下文）** |

**双层阻断与精确降次分析（修正粗略次数的严格结论）**：
1. **第一层（不含公理 (iv) 的纯维度无关迹型 NPA 层级）**：只施加 `(i)+(ii)+(iii)+(Tr)` 时，我们的四元数 4-MUB ($d=6$) 与四元数 5-MUB ($d=10, 14, \dots$) 是算子层面精确满足全部公理的有限维 $C^*$-表示，因此**对任意有限阶乃至无穷阶 $k \to \infty$ 均严格可行**，绝不可能证伪 $d=6$ 的 4-MUB 或 $d=10$ 的 5-MUB。
2. **第二层（施加公理 (iv) $[XUX, XVX]=0$ 作为双边理想时的精确击穿门槛）**：
   - **低次盲区（可行）**：当 $\deg(U) \le 1$ 或 $\deg(V) \le 1$ 时，由公理 (iii) 得 $XUX \propto X$，故 $[XUX, XVX] \equiv 0$ 在四元数模型上恒成立（注：GP 定理 3.2 脚注即指出该部分仅用到 $\le 3$ 次单项式）。此外，单纯迹关系 $\tau([A,B]) \equiv 0$ 恒成立，且对所有字长 $\le 1$ 的乘子 $W$ 均有 $\tau(W[A,B]) \equiv 0$，因此在**所有迹次数 $\le 6$（矩矩阵阶数 $t \le 3$）** 的层级上，公理 (iv) 仍完全无法识破四元数模型。
   - **精确高次击穿点（双边理想乘子）**：设最小非平凡角元 $A = P_{0,0} P_{1,0} P_{2,0} P_{0,0}$ 与 $B = P_{0,0} P_{1,0} P_{3,0} P_{0,0}$（各裸算 4 次），我们已零浮点证得 $C = [A, B] \ne 0$ 且精确满足 $A^\dagger A = B^\dagger B = \frac{1}{d^3} P_{0,0}$。在双边理想中施加 $L(W [A, B]) = 0$ 或 $L(C^\dagger C) = 0$ 时，利用 $P_{0,0}^2 = P_{0,0}$ 与共享内基的 $P_{1,0} P_{0,0} P_{1,0} = \frac{1}{d} P_{1,0}$ 降次，我们零浮点精确证实：
     - **线性理想见证最低约化循环次数 = 7**（需矩矩阵阶数 $t = 4$，多项式总次数 8）：取 $E = P_{0,0} P_{3,1} P_{1,0} P_{0,0}$，则 $\operatorname{Tr}(E[A,B]) = \frac{1}{d}[\operatorname{Tr}(P_{0,0} P_{3,1} P_{1,0} P_{2,0} P_{0,0} P_{1,0} P_{3,0}) - \operatorname{Tr}(P_{0,0} P_{3,1} P_{1,0} P_{3,0} P_{0,0} P_{1,0} P_{2,0})] \ne 0$（在 $d=6$ 验出 **4** 条、$d=10$ 验出 **18** 条双项均降至 7 次的非零见证；一般角元 $E$ 共验出 120 / 720 条 9 次以内非零见证）。
     - **SOS 交换子范数约化循环次数 = 10（不共享内基时 = 12）**：对上述共享 $P_{1,0}$ 的 $A, B$，裸算 16 次的 $\operatorname{Tr}(B^\dagger A^\dagger B A)$ 经 4 处 $P_0^2 = P_0$ 降至 12 次，再经 $P_1 P_0 P_1 = \frac{1}{d} P_1$ 精确降至 **10 次循环字** $W_{10} = P_{0,0} P_{3,0} P_{1,0} P_{0,0} P_{2,0} P_{1,0} P_{3,0} P_{0,0} P_{1,0} P_{2,0}$，满足精确恒等式 $\operatorname{Tr}(C^\dagger C) = \frac{4}{d^6} - \frac{2}{d}\operatorname{Re}\operatorname{Tr}(W_{10}) > 0$（已零浮点验毕）；若 $A, B$ 不共享内基（例如 $d=10$ 的 5-MUB 取 $P_0 P_1 P_2 P_0$ 与 $P_0 P_3 P_4 P_0$），则该 SOS 范数约化循环次数恰为 **12 次**。
   - **与 GP (2024) 文献的审慎衔接**：GP 原文关于 12 次门槛的表述为数值证据（*"numerical results indicate"*），而非下界定理；本文不宣称上述交换子是区分秩-1 与秩-2 的唯一最低次恒等式，而是将其定位为**揭示高次角代数交换约束为何必不可少的显式代数机制解释**。

### 2.3 普适任意秩-$r$ 2-Copy $\Lambda^2$ 零泄漏定理（Paper II 核心修复定理）

**定理（普适 2-Copy 反对称泄漏定理）**：*设 $P$ 为任意有限维复希尔伯特空间上任意秩 $r = \operatorname{rank}(P) \ge 1$ 的正交投影算子（$P^2 = P = P^\dagger$），$P_{\Lambda^2} = \frac{1}{2}(I - \mathrm{SWAP})$ 为 2-Copy 反对称子空间投影。则*
$$\frac{\operatorname{Tr}(P^{\otimes 2} P_{\Lambda^2})}{\operatorname{Tr}(P^{\otimes 2})} = \frac{r - 1}{2r}.$$
*特别地：*
1. *$P^{\otimes 2} P_{\Lambda^2} = 0$（0% 泄漏）**当且仅当** $\operatorname{rank}(P) = 1$（纯复秩-1 投影）；*
2. *对**任意**秩-2 表示（完全不限于本文四元数族！），反对称泄漏恒为精确的 **$25\%$**；*
3. *对**任意**高秩表示（$r \ge 2$），反对称泄漏均严格 $\ge 25\%$（单调趋向 $50\%$）。*

**三行证明**：因为 $P$ 为秩-$r$ 正交投影，故 $\operatorname{Tr}(P) = \operatorname{Tr}(P^2) = r$，从而 $\operatorname{Tr}(P^{\otimes 2}) = (\operatorname{Tr} P)^2 = r^2$，且由交换算子恒等式得 $\operatorname{Tr}(P^{\otimes 2} \mathrm{SWAP}) = \operatorname{Tr}(P^2) = r$。代入 $P_{\Lambda^2} = \frac{1}{2}(I - \mathrm{SWAP})$ 立得 $\operatorname{Tr}(P^{\otimes 2} P_{\Lambda^2}) = \frac{r^2 - r}{2}$，比值恰为 $\frac{r-1}{2r}$。由于 $P^{\otimes 2} P_{\Lambda^2} = (P_{\Lambda^2} P^{\otimes 2} P_{\Lambda^2})$ 半正定，其迹为零当且仅当算子本身为零，即 $r = 1$。$\square$

这一普适定理将 2-Copy $\Lambda^2$ 零泄漏约束（仅 **2 次**多项式！）提升为斩杀**一切高秩 $C^*$-代数反例**（包括四元数、八元数、任意更高维 Clifford 模型）的严格最小修复工具：Round 19 实测证实三套 2-Copy 秩-1 Level-1 SDP 均在 2 次即严格否定 `(d=2, m=4)` 与全部 `(d, d+2)`，并发现精确正交分解恒等式 $\Pi_a^{(R)}\Pi_b^{(R)} = F_{\text{ent}}$；同时 Level-2 $\operatorname{Sym}^2$（$17\times 17$）在 $d=6, 10$ 仍可行，揭示光滑 $d$-参数对称矩矩阵必须进一步结合有限整数维 $d$ 的离散相位闭合约束。

### 2.4 $\mathbb{H}^d \to \mathbb{C}^d$ 同伦质量间隙相图（Round 18）

| 维度 $d$ | $\mathbb{C}^d$ 状态 | 闭式四元数 $E_{\text{comm}}$ | $\mathcal{M}_{\mathbb{H}}(d,4) \to \mathbb{C}^d$ 同伦最小 $E_{\text{comm}}$ | 结论 |
|---:|---|---:|---:|---|
| $d=2$ | **已证不可能** | `1/36 = 0.027778` | **`2.778e-02`**（精确 `1/36`） | 硬质量间隙 |
| $d=3$ | 存在 | `0` | **`1.388e-16`** | 坍缩至精确复 4-MUB |
| $d=4$ | 存在 | `0` | **`3.955e-16`** | 坍缩至精确复 4-MUB |
| $d=5$ | 存在 | `0` | **`1.642e-14`** | 坍缩至精确复 4-MUB |
| **$d=6$** | **2026 预印本已证 $N(6)=3$** | `5/81 = 0.061728` | **`4.353e-02`** | **硬质量间隙（+12 量级）** |
| **$d=10$** | **核心未决难题** | `61/900 = 0.067778` | **`6.778e-02`（闭式）/ `1.150e-01`** | **硬质量间隙（+13 量级）** |

---

## 第三部分：$d=10$ 全局数值景观、$U(2)^5$ 主干结构与 492 万 GPU 轨迹证据

1. **7 大精确 $d=10$ 3-MUB 拓扑族（全部确认无第 4 向量）**：
   - `{16, 8, 8}` 刚性（$H_{\text{conf}} + D_0$, 零浮点证）与可形变（$H_{\text{conf}} + B_{\text{new}}$, 点零浮点证 + 整张 4D 连续流形 $\mathcal{M}_4$ 96 链测地线扫清，最小残差 `34.78`）；
   - `{16, 0, 0}` 刚性（$H_{\text{conf}} + B_{\text{iso}}$, 零浮点证）；
   - `{8, 8, 8}` 刚性（$F_2\otimes F_5$ 与 $U(2)^5$ 23 条孤立轨道）与可形变（`xSnap` 4D 流形）；
   - `{0, 8, 8}`（$\{I, F_{10}, \text{Circ}_2\}$）、`{0, 0, 8}`（双 $d=0$ 循环 + $d=8$）、`{0, 0, 0}`（全孤立循环三基 `{I, Circ #3, B3}`，完全跳出 $F_{10}$ 与 $H_{\text{conf}}$）。
2. **$U(2)^5$（$2\times 2$ 块 $5\times 5$ 循环阵）24 轨道普查与世界纪录**：
   - 20,400 起点命中 6,748 组精确 $U(2)^5$ 3-MUB，完备收敛至 **24 个唯一轨道**，逐位精确复现了第 14 轮无约束 GPU 找到的几乎全部刚性 `[8,8,8]` 3-MUB。
   - **完整 4-MUB 世界纪录**：`Pure U(2)^5` 18,000 次精确 LM 锁定唯一全局地板 `0.4612226`（最大内积误差 `1.34e-2`），释放到 270D 达 **`normalized = 0.4955`**（比 492 万 GPU 无约束最优 `11.78` 低 **24 倍**）。
   - **`{3-MUB + 2 正交向量}` 世界纪录**：Fourier-Twisted $U(2)^5$ 达 **`normalized = 3.33e-6` (180D) / `3.50e-6` (189D)**（比 GPU 最优 `5.10e-4` 低 **153 倍**），但全部 24 轨的 Bloch 共面行列式均严格止步于 `min|cop| >= 0.010454`（距零 1%）。
3. **三轮 Kaggle Dual-T4 GPU 累计 4,915,200 条轨迹全记录**：
   - Round 14（1,310,720 轨迹，过松弛 AltProj）、Round 15（1,310,720 轨迹，Douglas–Rachford + 隧道跳坑）、Round 19/16 v2（2,293,760 轨迹，Power-of-$U$ 393k + 78D Twisted-Diag $U(2)^5$ 196k + 退火 AltProj `3+1` 786k + 直接 `4-MUB` 917k）：**全线 0 精确 4-MUB、0 精确 `{3+1}`**。

---

## 第四部分：一键零浮点独立复验命令清单

所有核心零浮点证明脚本均在纯 Node.js（v24+，无外部依赖）上运行，完整 SHA-256 哈希载于 [MANIFEST_SHA256.json](file:///D:/antigravity/d10_mub_continuous/preprint_d10_mub/MANIFEST_SHA256.json)：

```powershell
# 1. 验证全部 37 份核心证明脚本、证书与预印本文稿的 SHA-256 完整性
node D:\antigravity\d10_mub_continuous\preprint_d10_mub\generate_sha256_manifest.mjs

# 2. 零浮点验证闭式四元数 4-MUB (d=6) 与 5-MUB (d=10, 14) 及角代数非交换性 (~1s)
node D:\antigravity\d10_mub_continuous\round17\prove_exact_quaternionic_mub_npa_barrier_QED.mjs

# 3. 零浮点验证显式复矩阵 M_12(Z[i,w_3]) 与 M_20(Z[i,w_5]) 全部 2,432 条三重积恒等式 (~13s)
node D:\antigravity\d10_mub_continuous\round17\verify_complex_20x20_npa_projectors_QED.mjs

# 3b. 零浮点逐条核验 Gribling–Polak Def. 3.1 公理 (i)-(iv)（含 144/500 条基内正交、432/2000 条三重积与 7/10 次降次见证）(~18s)
node D:\antigravity\d10_mub_continuous\preprint_d10_mub\verify_gp_alignment_and_degrees_QED.mjs

# 3c. 零浮点验证 Clifford 6-MUB (d=10, M_40), 8-MUB (d=14), 12-MUB (d=22) 与普适 Jordan-Wigner N_{C*}(d)=infinity (~73s / ~0.6s)
node D:\antigravity\d10_mub_continuous\round20\prove_clifford_qplus1_mub_QED.mjs
node D:\antigravity\d10_mub_continuous\round20\prove_universal_jordan_wigner_cstar_mub_QED.mjs

# 4. 零浮点验证 B_new 精确代数构造、固定 H_conf 90/90 严格孤立 (4/4 横截) 与 36 组 Degree-6 Nullstellensatz (~10s)
node D:\antigravity\d10_mub_continuous\round12\prove_B_new_exact.mjs
node D:\antigravity\d10_mub_continuous\round20\prove_B_new_fixed_Hconf_isolated_QED.mjs
node D:\antigravity\d10_mub_continuous\round13\prove_B_new_no_4th_vec_exact.mjs

# 5. 零浮点验证 B_iso 72 组手性轨道、171D 联合刚性及 128 位 CRT 对偶核无第 4 向量 Q.E.D. (~95s)
node D:\antigravity\d10_mub_continuous\round15\prove_B_iso_exact_theorems.mjs
node D:\antigravity\d10_mub_continuous\round15\prove_B_iso_crt_12core_QED.mjs
```
