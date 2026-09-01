# FLEKS 边界条件设计与迁移

本文是 FLEKS 边界条件（boundary condition，BC）体系的设计与重构总文档。

它取代 `PARAM.XML` 中 `#ABSORB` 条目引用的、实际上并不存在的
`Doc/Absorbing_BC_Implementation_Plan.md`。

相关源码：`include/BC.h`、`include/Particles.h`、`include/Pic.h`、
`src/Pic.cpp`、`src/Particles.cpp`、`src/Domain.cpp`。

---

## 1. 背景与动机

FLEKS 的边界条件原本由 `include/BC.h` 中一个 67 行的 `class BC` 承载：11 个
`static const int` 枚举值 + `amrex::IntVect lo, hi` + 一个字符串解析器
`num_type()`。**粒子边界（`#PARTICLEBOXBOUNDARY`）与电磁场边界
（`#FIELDBOXBOUNDARY`）共用同一个枚举和同一个解析器**，二者之间没有任何
"该类型是否属于本域"的校验。

这带来三类问题：

1. **类型错配被静默接受** —— 粒子侧可以写 `conducting`，场侧可以写 `reflect`，
   都能解析成功，但分派时没有任何分支处理它们。
2. **文档与代码不一致** —— `PARAM.XML` 把 `open` 列为合法选项，但
   `num_type()` 遇到 `open` 会直接 `Abort`。
3. **E 场与 B 场走不同路径** —— 所有 E 场调用都向 `apply_BC` 传
   `bc == nullptr`，导致 E 永远拿不到零梯度处理。

本文档给出缺陷清单、粒子域/场域的适用性对比、新的类型体系、输入命令表示、
向后兼容层与分阶段迁移路线。

---

## 2. 重构前的实现

### 2.1 数据结构

| 对象 | 位置 | 说明 |
| --- | --- | --- |
| `class BC` | `include/BC.h:17-64` | 枚举 + `lo`/`hi` + `num_type()` |
| `ParticlesInfo::pBCs` | `include/Particles.h:115` | `amrex::Vector<BC>(10)`，**按物种**索引，容量硬编码 10 |
| `Particles::bc` | `include/Particles.h:400` | 由 `set_bc()`（`:813`）从 `pBCs[iSpecies]` 注入 |
| `Pic::bcBField` | `include/Pic.h:260` | 全场**单一一份**，B 与 E 共用 |

> **`pBCs` 动态定尺寸的陷阱（实测踩坑）**
>
> 重构把 `ParticlesInfo::pBCs` 从固定容量 10 改为按 `#PLASMA nSpecies` 动态分配，
> 这暴露了一个隐藏的构造顺序依赖：**`Particles` 对象可能在 `pBCs` 定尺寸之前就被构造。**
>
> `src/TestParticles.cpp:9-13` 的 `make_test_particles_info()` 返回一个丢弃式的
> 默认 `ParticlesInfo`，其 `pBCs` **永远为空**；而 `Particles` 构造函数
> （`src/Particles.cpp:60`）会读 `pInfo.particle_bc(speciesID)`。固定容量为 10 时
> 这只是一次无害的默认读取，改成动态尺寸后变成越界访问，在
> `ParticleTracker::post_regrid()` 中直接段错误（SIGSEGV，测试报告 Errorcode 11）。
> **只有 `#PARTICLETRACKER T` 的测试会走到这条路径 —— 即 beam 与 photoionization。**
>
> 修复：`particle_bc()` 做越界保护，越界时返回全 `coupled` 的默认项
> （与相邻的 `initial_sup_id()` 已有的守卫风格一致）。
>
> **后续约定**：任何对 `pBCs` / `pBCsSet` 的新访问都必须同样加界检查。
> `Pic::post_process_param()` 里的 `resize(nSpecies)` 只保证 Pic 自己的 `pInfo`
> 被定尺寸，管不到 ParticleTracker 侧。

### 2.2 解析

`src/Pic.cpp:39-60`，两个命令都调用同一个 `BC::num_type()`。
命令路由在 `src/Domain.cpp:949-959`。

### 2.3 场侧分派

| 函数 | 位置 | 行为 |
| --- | --- | --- |
| `Pic::apply_BC` | `src/Pic.cpp:3595-3749` | 主入口；含 conducting > absorb > inflow 的**早退优先级链**；`bc == nullptr` 时走"全域 `func`"分支 |
| `Pic::use_float` | `include/Pic.h:637-687` | 零梯度拷贝；仅对 `outflow \|\| inflow \|\| fixed` 面生效 |
| `Pic::apply_conducting_wall` | `src/Pic.cpp:3752+` | PEC：B_n=0、B_t 镜像；E_t=0、E_n 镜像（`isB` 约定见 `:3802-3810`） |
| `Pic::apply_absorbing_wall` | `src/Pic.cpp:3820-3888` | 一阶匹配阻抗吸收；`dt <= 0` 时整体 return（`:3831-3833`） |
| `Pic::apply_inflow_wall` | `src/Pic.cpp:3891-3966` | 零梯度拷贝；`fixed` ≡ `inflow`；无 `#INFLOW` 块时 no-op |
| `Pic::apply_wave_field` | `src/Pic.cpp:4026-4096` | 硬源注入；`iField`：0=B，1=E |
| `Pic::apply_centerPlasma_BC` | `src/Pic.cpp:3970+` | 仅 hybrid；把离子矩镜像进墙外幽灵格 |

### 2.4 粒子侧分派

| 函数 | 位置 | 行为 |
| --- | --- | --- |
| `Particles::reflect_or_delete_particle` | `include/Particles.h:601-637` | 仅分派 `absorb`/`inflow`（删除 + 记账）与 `reflect`（镜面反射）；其余落到 `is_outside_active_region`；`iLev > 0` 时早退 |
| `Particles::absorb_tally` | `include/Particles.h:639-644` | 按面（2*d + {0=lo,1=hi}）累计 count/charge/mass |
| `Particles::inject_particles_at_boundary` | `src/Particles.cpp:~480-611` | 6 维布尔式枚举"不注入"的面类型，默认分支用 Maxwellian 重播种 |
| `Particles::inject_flux_at_inflow_faces` | `src/Particles.cpp:700-746` | 仅对 `inflow` 面；通量加权半空间 Maxwellian |
| 矩镜像折叠 | `src/Particles.cpp:1270-1320` | `doFold = reflect \|\| conducting \|\| inflow` |
| 电流（jHat）镜像折叠 | `src/Particles.cpp:2058-2117` | 同上判据 |
| 波速度扰动 | `src/Particles.cpp:343-366` | 读**粒子** BC 的 `wave` 标记 |

---

## 3. 缺陷清单

全部经实测代码确认。

| 编号 | 缺陷 | 位置 |
| --- | --- | --- |
| D1 | 粒子与场共用枚举，跨域错配被静默接受 | `include/BC.h:36-63` |
| D2 | `PARAM.XML` 文档化的 `open` 在 `num_type()` 中 Abort；`tests/` 下 0 个文件用 `open`（全用 `outflow`） | `PARAM.XML:1801,1834` vs `include/BC.h:59-62` |
| D3 | 粒子侧接受 `conducting`/`wave`/`fixed`，但 `reflect_or_delete_particle` 只分派 `absorb`/`inflow`/`reflect`，其余静默走"越界即删" | `include/Particles.h:601-637` |
| D4 | 场侧接受 `reflect`/`vacuum`，但 `apply_BC` 无对应分派，退化成 `func`（均匀态/流体接口填充） | `src/Pic.cpp:3595-3749` |
| D5 | **E/B 路径不对称**：所有 E 场调用传 `bc == nullptr`，永远拿不到 `use_float()` 零梯度；B 场走 `use_float` + `apply_inflow_wall` | `src/Pic.cpp:764, 2162-2164, 2775, 3250` vs `763/767` |
| D6 | `apply_BC` 是**早退优先级链**（conducting > absorb > inflow），混合面类型（x=conducting、y=absorb）在同一次调用中互相吞掉；实践中靠 `src/Pic.cpp:765-766`、`812-815` 手工连调绕过 | `src/Pic.cpp:3610-3637` |
| D7 | `BC::periodic` 与 `BC::coupled` **从不被分派**（全仓 0 命中）；写 `periodic` 但 `#PERIODICITY F` 时不报错、静默退化 | 全局 |
| D8 | `conducting` 在矩/电流折叠里是镜像面，在粒子删除路径里却不是反射面 | `src/Particles.cpp:1284-1286`、`2072-2074` vs `include/Particles.h:617-634` |
| D9 | 波粒子速度扰动判据读**粒子** BC，而 `tests/bc_wave/PARAM.in.*` 只设场边界，该路径不可达 | `src/Particles.cpp:343-366` |
| D10 | `PARAM.XML:1881` 引用 `Doc/Absorbing_BC_Implementation_Plan.md`，该文件不存在 | `PARAM.XML:1880-1881` |
| D11 | `src/Particles.cpp:558-566` 注释称 `"open"` 映射到 `BC::unset`，实际是 Abort，`:567-572` 为死代码 | `src/Particles.cpp:567-572` |
| D12 | 场命令参数名不统一：部分文件写 `fieldBoxBoundaryLo`，部分写 `fieldBoxBoundaryLo`。**仅影响可读性，不影响行为** —— `ReadParam::read_var()` 是**按位置**读取的，见 8.1 节 | `PARAM.XML:1832-1851` vs `src/Pic.cpp:56-57` |
| D13 | `#PERIODICITY` 文档参数名为 `isPeriodicX/Y/Z`，代码却在 `nDim` 循环里 `read_var("isPeriodic", …)` 三次 | `PARAM.XML:1283-1286` vs `src/Domain.cpp:1148-1153` |

---

## 4. 粒子域 vs 场域：适用性对比

| 输入字符串 | 粒子域 | 场域 | 物理含义与数值注意点 |
| --- | --- | --- | --- |
| `periodic` | 可用 | 可用 | 必须与 `#PERIODICITY` 一致；实际由 AMReX `FillBoundary(periodicity())` 处理。改为**自动派生**后成为可靠的按维早退标记（见第 5 节） |
| `coupled` | 可用（默认） | 可用（默认） | MHD-AEPIC 耦合：场取流体接口状态，粒子按 `fi` 重新播种。**从不显式分派**，是"未指定"的哨兵 |
| `outflow` / `open` | 可用 | 可用 | 粒子：跨越即删、不记账、不维持外部布居（Hybrid-VPIC "open"）。场：零梯度 Neumann 幽灵格（`use_float`）。`open` 为 `outflow` 的别名 |
| `vacuum` | 可用 | 可用（待实现） | 粒子：删除且不记账。场：外部无场。当前场侧无分派（D4），先登记为合法，实现留作后续 |
| `reflect` | 可用 | **不适用** | 粒子：镜面反射（位置 + 法向速度取反），能量守恒，适合对称面/静电壁。场侧请用 `conducting` 或 `symmetry` |
| `absorb` | 可用 | 可用 | **同名不同实现**：粒子侧 = 删除 + 按面记账（`absorb_tally`）；场侧 = 一阶匹配阻抗吸收（Engquist–Majda / Silver–Müller 离散形式，`decay`/`drive` 由 `#ABSORB charSpeed` 控制）。两者配对才是"完全吸收边界" |
| `inflow` | 可用 | 可用 | 粒子：通量加权半空间 Maxwellian 在物理面注入 + 删除出界粒子。场：上游 B + E 零梯度。**必须有 `#INFLOW` 块**，否则场侧静默退化为零梯度 |
| `thermal`（新增） | 可用 | 不适用 | 删除 + 按壁温半程 Maxwellian 再发射（等离子体-壁再循环）。当前未实现，先登记类型与解析 |
| `conducting` | **不适用** | 可用 | PEC 壁：B_n=0、B_t 镜像；E_t=0、E_n 镜像。**仅对静止导体正确** —— 流动等离子体需要 E_t = −u×B，错用会导致能量发散 |
| `symmetry`（新增） | 不适用 | 可用 | 镜像对称面：B_n=0 / B_t 镜像；E_n=0 / E_t 镜像。**与 PEC 的区别恰在 E_n**（PEC 镜像 E_n、置零 E_t），这是现有 `conducting` 无法表达的工况 |
| `fixed` | 不适用 | 可用（deprecated） | Dirichlet 钉住。**当前实现与 `inflow` 完全等价**（`use_float` 与 `apply_inflow_wall` 同等处理），标记 deprecated |
| `wave` | 不适用（仅作标记） | 可用 | 硬源/天线注入，需配套 `#WAVEBC`。粒子侧波速度扰动应由**场侧** `wave` 面驱动 |
| `unset` | 内部哨兵 | 内部哨兵 | 不对用户开放 |

### 4.1 错配矩阵

| 写法 | 处理 |
| --- | --- |
| 粒子侧写 `conducting` / `wave` / `fixed` / `symmetry` | 弃用警告 + 按别名表映射 |
| 场侧写 `reflect` / `thermal` | 弃用警告 + 按别名表映射 |
| 任一侧写未知字符串 | **始终 Abort**，错误消息列出该域全部合法值 |

---

## 5. 确认项一：`#PERIODICITY` 是否必要？

**结论：不可移除。但边界命令里的 `periodic` 取值目前 100% 惰性，应与几何
自动对齐。**

### 5.1 为什么不能移除

**(a) 读取顺序是硬约束。** 参数是两遍读的：

- `Domain::Domain()` → `prepare_grid_info()`（`src/Domain.cpp:38`）
  → `read_param(true)`（`src/Domain.cpp:262`）。`#PERIODICITY` 属于 grid 命令
  （列表见 `src/Domain.cpp:913-918`）。
- 紧接着 `gm.define(centerBox, &domainRange, coord, periodicity.getVect())`
  （`src/Domain.cpp:303`）—— AMReX `Geometry` 在此拿到周期性。
- 场边界命令要到 `src/Domain.cpp:83` 的 `read_param(false)` 才被解析，而
  `pic` 在 `:70-75` 才构造。

**即：Geometry 必须在场边界字符串被读到之前就已确定周期性。**

**(b) 它是拓扑属性，不是边界条件。** `Geom(iLev).isPeriodic(d)` 被约 70 处消费：

- `FillBoundary` / `SumBoundary` 的幽灵交换与粒子矩折叠
  （`src/Particles.cpp:1338`、`1514`）
- `apply_float_boundary` 的循环边界（`src/GridUtility.cpp:255-285`，非周期维
  要多扫一行）
- MHD→PIC 取点时的 `shift_periodic_index`（`src/FluidInterface.cpp:606-608`）
- 绘图几何重建（`src/PicIO.cpp:1043-1047`）与周期节点裁剪
  （`src/PicIO.cpp:255-260`、`318-326`）
- AMReX `BoxArray` / `DistributionMapping` 构造

**(c) `isFake2D` 强制覆盖。** `src/Domain.cpp:265-267`：

```
// If MHD is 2D, PIC has to be periodic in the z-direction.
if (isFake2D)
  set_periodicity(iz_, true);
```

与任何边界设置无关。

### 5.2 为什么边界命令里的 `periodic` 是惰性的

`Grid::update_cell_status`（`src/Grid.cpp:236-305`）：

1. fabbox 全部单元先标 `bit::set_lev_boundary`（`:244-252`）；
2. valid 单元被清成 `set_not_lev_boundary`（`:254-270`）；
3. `cellStatus[iLev].FillBoundary(Geom(iLev).periodicity())`（`:290`）
   —— **在周期方向把对面 valid 单元的 `not_lev_boundary` 标志搬进幽灵环**；
4. 剩余 `is_lev_boundary` 单元再做 `is_inside_domain` 判定（`:293-305`）。

因此边界算子的守卫 `if (!bit::is_lev_boundary(...)) return;` 在周期方向的幽灵
格上恒为假，`BC::periodic` 永远不会被分派（全仓 0 命中）。

### 5.3 采纳方案：自动派生

`#PERIODICITY` 保留为拓扑的唯一权威。Geometry 确定之后，凡
`Geom.isPeriodic(d)` 的维，强制

```
bcField.lo[d] = bcField.hi[d]   = FieldBC::periodic
pBCs[i].lo[d] = pBCs[i].hi[d]   = ParticleBC::periodic   （各物种）
```

用户写了冲突类型则告警。收益：

- 存储值不再能与 Geometry 矛盾（修 D7）；
- `apply_field_bc` 可以**按维**早退，而不必依赖 `isAllPeriodic()` 这个粗粒度
  特例。

**落点约束**：必须在 `gm.define()`（`src/Domain.cpp:303`）之后、
`read_param(false)`（`:83`）消费 BC 之前执行。

### 5.4 已否决的方案

- **方案 B**（把场边界命令移入 grid 遍，用 Domain 级容器暂存原始字符串，
  第二遍再解析）—— 侵入大，且绕不开 5.1 的 (b)(c) 两点，`#PERIODICITY` 仍需存在。
- **方案 C**（仅做一致性校验）—— 冗余依旧保留，用户仍要在两处重复声明。

---

## 6. 确认项二：全 PIC 与 hybrid 的场边界一致性

### 6.1 实测差异

|  | 全 PIC（`solveEM=T`） | hybrid（`useHybridPIC=T`，`src/Pic.cpp:286-287` 强制 `solveEM=false`） |
| --- | --- | --- |
| 演化场 | `nodeE`（隐式 GMRES，`update_E_impl` `src/Pic.cpp:2127`）+ `nodeB`（Faraday） | **仅 `centerB`**（RK4/SSPRK3 Faraday，`update_B_hybrid` `src/Pic.cpp:2914`） |
| E 的地位 | **独立状态量** | **导出量**：`assemble_ohm_E(centerBin, centerBtimeAvg, Eout, iLev, hstep)`（`src/Pic.cpp:2602`） |
| E 的 BC 调用点 | `nodeE`/`nodeEth`：`apply_BC(..., &Pic::get_node_E, iLev)` **bc=nullptr**（`:2162-2164`）+ conducting（`:2165-2168`）+ absorbing（`:2169-2172`）+ wave（`:2173-2179`）。**缺 `apply_inflow_wall`** | `Eout`：`apply_BC(..., &Pic::get_center_E, iLev)`（`:2775`）+ conducting（`:2776`）+ absorbing（`:2778`）+ **`apply_inflow_wall`**（`:2782`）；`centerEhybrid`：`FillBoundary` + `apply_BC`（`:3249-3251`）+ **`apply_inflow_wall`**（`:3256-3257`） |
| `nodeE` 角色 | 状态量 | 仅 IC/restart 种子：`average_node_to_center(nodeE[iLev], centerEhybrid[iLev])`（`src/Pic.cpp:807`） |

即：hybrid 里 E 是 B 与矩的泛函，E 的"边界条件"实质是**幽灵环闭包**（让边格
Faraday 旋度与 Boris push gather 有定义），不是给独立自由度加约束；全 PIC 里
E 是真状态量。**且当前全 PIC 的 `nodeE` 恰好缺了 hybrid 拥有的
`apply_inflow_wall`。**

### 6.2 统一设计：一个物理面一份类型，两个应用算子

用户只给出**一份**面类型，`isB` 选择 B 算子 / E 算子。

支持该简化的实测证据 —— "一份类型、两个算子"已是既成事实：

```2761:2768:src/Pic.cpp
    apply_BC(cellStatus[iLev], centerLapB[iLev], 0, centerLapB[iLev].nComp(),
             &Pic::get_center_B, iLev, &bcBField);

    curl_center_to_center(centerLapB[iLev], centerHyperE[iLev],
                          Geom(iLev).InvCellSize());
    centerHyperE[iLev].FillBoundary(Geom(iLev).periodicity());
    apply_BC(cellStatus[iLev], centerHyperE[iLev], 0,
             centerHyperE[iLev].nComp(), &Pic::get_center_E, iLev, &bcBField);
```

同一个 `bcBField` 对象既喂 B 又喂 E。而且所有墙算子都已内置 E/B 区分参数：
`apply_conducting_wall(..., bool isB)`、`apply_absorbing_wall(..., bool isB)`、
`apply_inflow_wall(..., bool isB)`、`apply_wave_field(..., int iField /*0=B,1=E*/, ...)`。
**E 与 B 的差异已经编码在算子内部，不需要两套输入。**

连带确认：D5 的根因是 `apply_BC` 对 E 传 `bc == nullptr`，而非缺少独立的
E 边界类型；单命令方案同样能修。

### 6.3 语义表

| 求解器 | `#FIELDBOXBOUNDARY` 的含义 |
| --- | --- |
| 全 PIC（`solveEM=T`） | 该类型同时约束演化量 E 与 B；B 算子和 E 算子按 `isB` 各取约定 |
| hybrid（`useHybridPIC=T`） | 该类型约束唯一演化场 `centerB`；E 由 `assemble_ohm_E` 导出，同一类型仅用于 E 幽灵环闭包 |

选项集合**不按求解器切分** —— 所有类型在两个求解器下都有对应算子落点。

---

## 7. 新类型体系

### 7.1 命名规范

| 层面 | 规范 | 示例 |
| --- | --- | --- |
| C++ 命名空间 | `ParticleBC` / `FieldBC`，大驼峰 | `FieldBC::conducting` |
| C++ 枚举成员 | 小驼峰 | `ParticleBC::thermal` |
| 输入字符串 | 全小写单词（保留历史拼写） | `outflow`、`conducting` |
| 存储容器 | `BoxBC<EnumT>` 模板 | `BoxBC<FieldBC::Type>` |
| 成员变量 | 小驼峰 | `bcField`、`pBCs` |
| 错误消息 | 带域名前缀，列出该域全部合法值 | `unrecognized particle boundary type 'conducting' …` |

### 7.2 数据结构

```cpp
// include/BC.h —— 统一的 lo/hi 存储容器，替换原 class BC
template <typename EnumT> struct BoxBC {
  amrex::IntVect lo, hi;
  BoxBC() : lo(EnumT::coupled), hi(EnumT::coupled) {}
  int face(int d, int side) const { return side == 0 ? lo[d] : hi[d]; }
};

// 解析契约：strict=true 时 legacy 别名直接 Abort；
// 否则返回映射值，并通过 warn 输出一次性弃用提示
namespace ParticleBC {
  enum Type { unset=-1, periodic=0, coupled=1, outflow=2, vacuum=3,
              reflect=4, absorb=5, inflow=6, thermal=7 };
  Type parse(const std::string& s, bool strict, std::string& warn);
}
namespace FieldBC {
  enum Type { unset=-1, periodic=0, coupled=1, outflow=2, vacuum=3,
              conducting=4, symmetry=5, absorb=6, inflow=7, fixed=8, wave=9 };
  Type parse(const std::string& s, bool strict, std::string& warn);
}
```

**为什么用 `namespace + unscoped enum` 而不是 `enum class`**：现有代码把类型值
直接存进 `amrex::IntVect`，并在约 60 处写作 `bc.lo[d] == BC::outflow`。非限定
枚举保留到 `int` 的隐式转换，使存储层与比较层零改动，把 diff 收敛到类型声明、
解析表与分派判据三处。`enum class` 需要 `static_cast` 全量改写，且无法直接
存入 `IntVect`。

枚举数值重新编号（不与旧值对齐）：这些值从不持久化到输入文件或重启文件，
重排没有兼容代价。

### 7.3 统一分派入口

```cpp
void Pic::apply_field_bc(const amrex::iMultiFab& status, amrex::MultiFab& mf,
                         int iStart, int nComp, GETVALUE func, int iLev,
                         bool isB);
```

取代现有样板：

```cpp
// 现状（src/Pic.cpp:2162-2172，全 PIC nodeE）
apply_BC(nodeStatus[iLev], nodeE[iLev], 0, nDim3, &Pic::get_node_E, iLev);  // bc == nullptr
apply_conducting_wall(nodeStatus[iLev], nodeE[iLev], 0, nDim3, iLev, bcBField, false);
apply_absorbing_wall(nodeStatus[iLev], nodeE[iLev], 0, nDim3, iLev, bcBField, false);

// 目标
apply_field_bc(nodeStatus[iLev], nodeE[iLev], 0, nDim3,
               &Pic::get_node_E, iLev, /*isB=*/false);
```

`apply_field_bc` 内部：

1. 始终使用 `bcField`（单一对象），按 `isB` 选择算子语义；
2. 把 `src/Pic.cpp:3610-3637` 的**早退优先级链换成逐面 switch**，混合面类型
   一次扫完（修 D6）；
3. 保留各墙算子内部的 `isB` 约定；`apply_wave_field` 的 `iField` 由
   `isB ? 0 : 1` 得到；
4. `bc == nullptr` hack 消失（修 D5）；
5. 全 PIC `nodeE` 自动获得 hybrid 已有的 `apply_inflow_wall`。

两个求解器读同一对象：

- 全 PIC：`apply_field_bc(nodeStatus, nodeB, ..., isB=true)` /
  `apply_field_bc(nodeStatus, nodeE/nodeEth, ..., isB=false)`
- hybrid：`apply_field_bc(cellStatus, centerB, ..., isB=true)` 作用于演化场；
  `apply_field_bc(cellStatus, centerEhybrid/centerEstage, ..., isB=false)`
  紧跟 `assemble_ohm_E`，语义注明为"幽灵环闭包，非独立约束"

---

## 8. 输入命令表示

```
#PERIODICITY            [保留，拓扑唯一权威 — 见第 5 节]
isPeriodicX/Y/Z       logical, default F

#PARTICLEBOXBOUNDARY    (per species, multiple=T)  [扩展]
particleBoxBoundaryLo/Hi   select: periodic | coupled | outflow | open | vacuum
                                   | reflect | absorb | inflow | thermal

#FIELDBOXBOUNDARY       [新增，场边界的唯一入口]
fieldBoxBoundaryLo/Hi      select: periodic | coupled | outflow | open | vacuum
                                   | conducting | symmetry | absorb
                                   | inflow | fixed | wave

#FIELDBOXBOUNDARY      [deprecated，#FIELDBOXBOUNDARY 的别名]

#BCSTRICT               [新增，可选]
strictBC                   logical, default F
                           T 时把所有弃用警告与一致性警告升级为 Abort

#ABSORB / #INFLOW / #WAVEBC   [保持不变，仅更新归属文档]
```

要点：

- **参数名与粒子侧对称**：`#PARTICLEBOXBOUNDARY` → `particleBoxBoundaryLo/Hi`；
  `#FIELDBOXBOUNDARY` → `fieldBoxBoundaryLo/Hi`。顺带修 D12。

### 8.1 参数名是位置性的，无需别名层

`ReadParam::read_var()` **按位置**读取，不做名称匹配 —— `description` 参数
只用于回显（`share/Library/src/ReadParam.h:153-174` 的模板版本直接
`ss >> var`；`:201-245` 的字符串版本读整行后取首词）。

推论：

- `fieldBoxBoundaryLo` / `fieldBoxBoundaryLo` / `fieldBoxBoundaryLo` 三者完全等
  价，旧输入文件无需改动即可继续工作。
- 重构时只需保证**读取的行数不变**（每个方向 lo/hi 各一行），参数名可以自
  由改名。
- 只有 `ReadParam::read_optional()`（`:258-323`）才真正比对名称，因此可选
  参数不能这样随意改名。

这是 D12 只需改文档、不需改代码的原因。
- **`#FIELDBOXBOUNDARY` 走代码层 deprecated 别名**（`PARAM.XML` 的 `alias`
  属性仅用于 `FLEKS0/1/2` 网格后缀，不承载 legacy 同义词）。
- 命令路由在 `src/Domain.cpp:949-959` 的 pic 命令列表中追加 `#FIELDBOXBOUNDARY`
  与 `#BCSTRICT`，`#FIELDBOXBOUNDARY` 保留。

---

## 9. 向后兼容层

| 上下文 | legacy 串 / 名 | 映射为 | 处理 |
| --- | --- | --- | --- |
| 命令名 | `#FIELDBOXBOUNDARY` | `#FIELDBOXBOUNDARY` | 弃用警告 |
| 参数名 | `fieldBoxBoundaryLo/Hi`、`fieldBoxBoundaryLo/Hi` | `fieldBoxBoundaryLo/Hi` | 静默接受（修 D12） |
| 两侧 | `open` | `outflow` | 静默接受（PARAM.XML 已文档化，属修 D2） |
| 两侧 | `vacume` | `vacuum` | 静默（保留历史拼写） |
| 粒子 | `conducting` | `reflect` | 弃用警告（与 `tests/bc_field` 现有手写法一致） |
| 粒子 | `fixed` | `inflow` | 弃用警告 |
| 粒子 | `wave` | `outflow` + 标记该面为波驱动面 | 弃用警告 |
| 场 | `reflect` | `conducting` | 弃用警告 |
| 场 | `fixed` | `inflow`（保持当前等价行为） | 弃用警告，后续拆为独立 Dirichlet |
| 两侧 | `unset` / 未知串 | — | 始终 Abort，错误消息列出该域全部合法值 |

弃用警告经 `amrex::Print()` 输出，每条仅打印一次（按"面 + legacy 串"去重）。
`#BCSTRICT T` 时全部升级为 `amrex::Abort`。

---

## 10. 跨域一致性校验 `validate_bc_pairing()`

在参数读取完成后（`Domain::read_param` 末尾）调用，逐面检查，默认告警、
`#BCSTRICT` 下 Abort：

1. `#PERIODICITY` 该维为 `T` 但边界写了非 `periodic` 类型
2. 边界写 `periodic` 但该维 `#PERIODICITY` 为 `F`
3. 同一维 lo/hi 一侧 `periodic` 另一侧非 `periodic`
4. 粒子 `inflow` 但场非 `inflow`/`fixed`
5. 场 `inflow` 但无 `#INFLOW` 块（当前静默退化）
6. 场 `wave` 但 `waveBC` 未激活（无 `#WAVEBC`）
7. 场 `conducting`/`symmetry` 但粒子为 `outflow`/`vacuum`/`absorb`（粒子穿墙）
8. 粒子 `absorb` 但场非 `absorb`（仅半吸收，提示）
9. hybrid 且场为 `conducting`/`symmetry` → 提示 E 侧仅为幽灵环闭包、
   `inflow` 需 `#INFLOW`

---

## 11. 迁移路线（五个阶段）

- **阶段 1**：拆分枚举 + 兼容层 + `open` 修复 + 周期性自动派生 + 一致性校验。
  **所有现有 `PARAM.in` 不改一行即可运行，数值结果逐位不变。**
- **阶段 2**：更新 `PARAM.XML`（新增 `#FIELDBOXBOUNDARY` / `#BCSTRICT`，
  标注 `#FIELDBOXBOUNDARY` deprecated，修正文档指向，修 D13 参数名）。
- **阶段 3**：机械迁移 `tests/` 下 15+ 个 `PARAM.in` 到 `#FIELDBOXBOUNDARY` +
  `fieldBoxBoundaryLo/Hi`。
- **阶段 4**：测试侧打开 `#BCSTRICT T`，确认零警告后固化。
- **阶段 5（后续独立任务）**：实装 `apply_field_bc` 统一分派（修 D5/D6）、
  新增 `apply_symmetry_wall`、实现 `vacuum`/`thermal` 的场/粒子侧行为、
  删除 legacy 别名表。

**阶段 1-4 保证行为不变；阶段 5 才触碰数值行为，需独立评审与回归。**

---

## 12. 附录：引用点清单

重构前的精确统计（`include/` + `src/`；`srcInterface/` 零引用）。

| 模式 | 文件 | 行数 | 域 |
| --- | --- | --- | --- |
| `bcBField` | `src/Pic.cpp` | 38 | 场 |
| `bcBField` | `include/Pic.h` | 1（`:260` 声明） | 场 |
| `BC::<member>` | `include/Pic.h` | 13 | 场 |
| `BC::<member>` | `src/Pic.cpp` | 10 | 场 |
| `BC::<member>` | `include/Particles.h` | 5 | 粒子 |
| `BC::<member>` | `src/Particles.cpp` | 14 | 粒子 |
| `bc.<member>`（实例限定） | `src/Particles.cpp:517-580` | 36 | 粒子 |

**注意**：`src/Particles.cpp:517-580` 的 36 处使用**实例限定**访问
（`bc.inflow`、`bc.outflow`、`bc.vacuum`、`bc.absorb`、`bc.unset`、`bc.reflect`），
而非 `BC::inflow`。它们与 `BC::` 形式一样需要一并改写。

**同名陷阱**：`src/Pic.cpp` 中的形参名 `bc` 是**场**侧（`bcBField` 传入），
`src/Particles.cpp` / `include/Particles.h` 中的成员 `bc` 是**粒子**侧。
批量替换时必须按文件区分，不可套用同一规则。
