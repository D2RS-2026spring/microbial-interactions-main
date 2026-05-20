📊 海洋微生物群落演替与互作网络动态特征复现项目
本项目基于 Quarto Book 技术栈构建，完整复现了原论文中涉及群落结构时间序列、网络整体特征（促进 vs 竞争）以及环境驱动力分析的三大核心生态学结论：

图 1A：微生物群落结构随时间演替的长期时间序列变化。

图 2A-C：微生物互作网络整体特征与动态变化趋势。

图 4：核心微生物互作网络随外部水温演变的环境依赖性综合面板（A-E）。

⚙️ 1. 环境准备 (Prerequisites)
本项目使用 R 语言 进行数据流式处理与绘图。请确保您的系统已安装 Quarto CLI 以及 R 语言运行环境。

1.1 安装基础工具
Quarto CLI: v1.3 或更高版本（官方下载链接）

R 语言: v4.0 或更高版本

1.2 安装 R 依赖包
为了保障代码在不同电脑上都能无缝运行，我们在代码中内置了清华大学 CRAN 镜像源，并强制采用预编译的二进制（binary）版本进行安装，完美绕过本地本地 C++ 静态编译引发的各类环境报错。

打开 R 终端或 RStudio，直接运行以下命令安装项目所需的全部核心扩展包：

R
install.packages(c("tidyverse", "data.table", "ggplot2", "ggpubr", "ggridges", "ggtext", "ggpattern", "moments"))
📁 2. 严谨的目录结构核对 (Directory)
本项目的所有分析脚本与数据读取均采用相对路径。请在运行前严格核对您的本地目录是否与下方完全一致：

Plaintext
MICROBIAL-INTERACTION-xxx/ (本地克隆的仓库根目录)
│
├── README.md                    # 👈 当前说明文档（请存放在这里）
├── LICENSE                      # 许可证文件
├── r_environment.Rproj          # RStudio 项目环境
├── .gitignore                   # Git 忽略文件（已自动忽略 .quarto 与 _book 缓存）
│
├── data/                        # 📜 基础元数据与群落丰度数据
│   ├── taxa_information.csv
│   └── data_sequences_0.1_rel_ab_0.5_occ_binned_4_days_with_temperature.csv
│
├── model_out/                   # 💾 网络模型输出（包含 1GB 巨型数据文件）
│   ├── MDR_smap_interactions_Scripps_Pier_ASVs_tp_2_28092024.csv
│   ├── strongest_interactions_Scripps_Pier_ASVs_tp_2_28092024.csv
│   └── keystone_taxa_Scripps_Pier_ASVs_tp_2_28092024.csv
│
├── scripts/                     # ⚙️ 原作者提供的原始 R 脚本
│   └── 4-Fig_2A-C_microbial_interactions.R
│
├── plots/                       # 🖼️ 核心科研图表导出目标文件夹
│
└── LL/                          # 📂 Quarto 项目核心工作区
    ├── _quarto.yml              # Quarto 配置文件
    ├── index.qmd                # 本地电子书首页说明
    ├── intro.qmd                # 引言说明
    ├── summary.qmd              # 总结说明
    ├── references.qmd / .bib    # 文献引用配置
    ├── 01-community-timeseries.qmd # 【重构分析】图 1A 与 图 4 环境依赖性代码
    └── 02-network-interactions.qmd # 【沙盒调用】图 2A-C 微生物互作整体网络特征代码
⚠️ 运行必看：请确保项目根目录下存在空的 plots/ 文件夹。代码中所有读取数据的相对路径均以 LL/ 文件夹为基准向上查找（如 ../data/、../model_out/），切勿随意更改外层文件夹的名称！

🚀 3. 一键复现与发布步骤 (Execution & Publish)
打开您的终端（Terminal 或 CMD），进入到 LL 子文件夹中进行编译：

📌 步骤 A：本地一键编译整个项目 (Render)
Bash
# 1. 必须先进入到包含 _quarto.yml 的 LL 工作区
cd LL

# 2. 一键编译整本 Quarto Book
quarto render
运行机制：该命令会流式调用 01-community-timeseries.qmd（自主重构计算）和 02-network-interactions.qmd（通过 local = TRUE 沙盒安全调用外层 ../scripts/ 下的原作者官方绘图流水线）。

📌 步骤 B：一键发布到云端 (Publish)
如果您需要将复现生成的电子书网页部署到 GitHub Pages 上供在线浏览，请在 LL 目录下直接运行：

Bash
quarto publish gh-pages
确认提示后输入 Y 即可。发布成功后，终端末尾会生成一个在线公网链接，老师点击即可直接查阅。

⏳ 硬件性能与耗时说明：
由于系统在编译时，两个 .qmd 文件均需要高频、流式地从外层 ../model_out/ 中读取磁盘上 1GB 左右的巨型网络系数 CSV 文件。根据您电脑固态硬盘（SSD）的读写性能，代码在执行到 fread 区块时，可能会卡顿 10 ~ 30 秒，这属于正常的数据 I/O 现象，系统并未死机，请耐心等待其全部编译完成。

🎯 4. 预期复现成果 (Expected Output)
4.1 本地/在线静态网页
编译完成后，会在 LL/ 目录下生成 _book/ 文件夹。其中包含：

01-community-timeseries.html：直观展示群落结构时间序列堆叠图、水温波动折线、以及图 4 综合面板。

02-network-interactions.html：展示群落“相亲相爱”与“互相内卷”整体网络特征的中文生态学解读。

4.2 出版级高清矢量图文件
在项目外层的 plots/ 文件夹下，会自动安全导出以下高清科研 PDF 图表：

Fig_1A_microbial_community_time_series.pdf（图 1A）

Fig_4_temperature_dependency.pdf（图 4 综合大图）

Fig_2A.pdf、Fig_2B.pdf、Fig_2C.pdf（由作者原生脚本生成的图 2 核心面板）
