# RS 构型学习工具

> 面向 315 化学（农）考研的手性分子 R/S 构型判断 3D 可视化学习工具

[![在线体验](https://img.shields.io/badge/在线体验-ad429w.kicp.fun-brightgreen)](https://ad429w.kicp.fun/)
[![3Dmol.js](https://img.shields.io/badge/3Dmol.js-2.0-blue)](https://3dmol.org/)
[![License](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

## 项目简介

考研有机化学中，手性分子的 R/S 构型判断是高频考点，也是难点。很多同学面对试卷上的楔形键结构式时，难以在脑海中建立三维空间关系，导致判断出错。

本工具通过 **3D 交互式可视化**，将平面结构式还原为可旋转、可缩放的三维分子模型，帮助理解基团的空间排布和 CIP 优先级排序，从而掌握 R/S 构型的判断方法。

工具同时提供**练习模式**，随机生成化合物进行构型判断训练，答错后显示详细解析（含 CIP 排序过程和判断步骤）。

## 功能特性

### 学习模式

| 功能 | 说明 |
|------|------|
| 预设化合物库 | 收录 30+ 考研高频化合物，覆盖乳酸、甘油醛、氨基酸、卤代烃、醇、环状化合物等 |
| 二维构建器 | 遇到没名字的题目？手动选择基团和键型（纸面/实楔朝前/虚楔朝后），还原试卷结构后生成 3D 模型 |
| SMILES 输入 | 支持通过 SMILES 字符串快速构建分子 |
| 单手性/双手性中心 | 支持含一个或两个手性中心的分子，覆盖酒石酸等典型考题 |
| 3D 交互 | 鼠标拖拽旋转、滚轮缩放，可切换氢原子显示、基团标签、构型标记 |
| 对映体切换 | 一键生成对映异构体，对比 R/S 空间差异 |
| Fischer/Newman 投影 | 支持查看 Fischer 投影式和 Newman 投影式 |

### 练习模式

| 功能 | 说明 |
|------|------|
| 随机出题 | 随机生成化合物 3D 模型，判断手性中心构型 |
| 多视角辅助 | 可旋转模型、切换判断视角辅助观察 |
| 详细解析 | 答错后显示 CIP 优先级排序和完整判断过程 |
| 答题统计 | 记录答题正确率，针对性强化 |

## 技术栈

| 类别 | 技术 | 说明 |
|------|------|------|
| 3D 渲染 | [3Dmol.js](https://3dmol.org/) | 分子三维可视化引擎 |
| 前端 | 原生 HTML/CSS/JavaScript | 单文件应用，无构建依赖 |
| 样式 | CSS 变量 + 响应式布局 | 深色主题，适配桌面和移动端 |
| 算法 | CIP 优先级排序 | BFS 逐层比较原子序数，支持双键/三键/苯环 |
| 验证 | Python | 科学正确性验证脚本 |

## 快速开始

### 方式一：在线使用

直接访问：**https://ad429w.kicp.fun/**

### 方式二：本地运行

本项目是纯前端单文件应用，无需安装任何依赖：

```bash
# 克隆项目
git clone https://github.com/yuyann-dev/rs-configuration-learning-tool.git

# 直接用浏览器打开
cd rs-configuration-learning-tool
# 双击 index.html，或用任意静态服务器打开
```

也可以用 Python 启动一个简单的本地服务器：

```bash
python -m http.server 8080
# 然后访问 http://localhost:8080
```

## 使用说明

### 学习模式

1. **选择化合物**：左侧面板从预设化合物库中点击选择，或使用「二维构建器」手动构建
2. **观察 3D 模型**：中间区域拖拽旋转、滚轮缩放
3. **判断构型**：根据 CIP 规则，将最小基团置于远离观察者的位置，观察其余三个基团的排列方向
4. **查看解析**：右侧面板显示 CIP 排序结果和构型标记

### 练习模式

1. 点击顶部「练习模式」切换
2. 系统随机生成化合物，选择 R 或 S
3. 答错后查看详细解析，理解判断过程
4. 可旋转模型从不同角度观察

### 二维构建器使用

1. 选择手性中心数量（1个或2个）
2. 为每个位置选择基团（下拉菜单）
3. 选择键型：纸面内 / 实楔朝前 / 虚楔朝后
4. 点击「生成 3D 模型」
5. 可点击「对映体」一键切换构型

## 目录结构

```
rs-configuration-learning-tool/
├── index.html              # 主应用（单文件，含HTML/CSS/JS）
├── scripts/                # 科学验证脚本
│   ├── verify_rs.py            # CIP排序与RS构型验证
│   ├── verify_rs_round2.py     # 第二轮验证
│   ├── verify_rs_round3.py     # 第三轮验证
│   ├── verify_round1.py        # 基础题验证
│   ├── verify_round2.py        # 进阶题验证
│   └── verify_exam_questions.py # 真题风格题验证
├── screenshots/            # 截图
│   ├── desktop.jpg             # 桌面端界面
│   └── mobile.jpg              # 移动端界面
└── README.md
```

## 科学验证

项目包含 Python 验证脚本，用于确保 CIP 优先级排序和 R/S 构型判断的科学正确性：

```bash
cd scripts
python verify_rs.py
python verify_exam_questions.py
```

验证覆盖：
- 基础手性分子（乳酸、甘油醛、2-氯丁烷等）
- 氨基酸系列（丙氨酸、丝氨酸、半胱氨酸、苯丙氨酸等）
- 双手性中心分子（酒石酸、苏阿糖等）
- 环状化合物（环丙烷衍生物、环己酮衍生物）
- 315 化学（农）真题风格题目

## CIP 优先级规则实现

工具内置完整的 Cahn-Ingold-Prelog 优先级排序算法：

1. **直接比较**：比较与手性中心直接相连的原子序数
2. **逐层 BFS**：原子序数相同时，向外层扩展，逐层比较
3. **重键处理**：双键视为连接两个相同原子，三键视为三个
4. **苯环处理**：芳环碳按 Kekulé 结构处理
5. **同位素**：D > H（质量数优先）

排序完成后，将最小基团置于远离观察者位置，其余三基团顺时针为 R，逆时针为 S。

## 部署

### 静态托管

本项目为纯静态应用，可部署到任意静态托管平台：

- **Vercel / Netlify**：直接连接 GitHub 仓库，自动部署
- **GitHub Pages**：推送到 `gh-pages` 分支
- **Nginx**：将 `index.html` 放入网站根目录

### Nginx 配置示例

```nginx
server {
    listen 80;
    server_name your-domain.com;
    root /var/www/rs-configuration-learning-tool;
    index index.html;

    location / {
        try_files $uri $uri/ /index.html;
    }
}
```

## 适用场景

- 315 化学（农）考研备考
- 有机化学手性分子教学
- 化学相关专业学生自学
- 化学竞赛训练

## 许可证

MIT License
