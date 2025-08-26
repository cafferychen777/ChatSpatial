#!/bin/bash

# Claude Code LaTeX Helper Commands
# 为 Claude Code 提供的 LaTeX 辅助命令集合

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
MAGENTA='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m'

# Claude Code 专用提示词生成器
generate_latex_prompt() {
    local task=$1
    local context=$2
    
    case $task in
        "formula")
            echo "Help me write a LaTeX formula for: $context"
            echo "Please provide the formula with proper LaTeX syntax and explain any special packages needed."
            ;;
        "table")
            echo "Create a LaTeX table with the following specifications: $context"
            echo "Include proper formatting, alignment, and any necessary packages."
            ;;
        "figure")
            echo "Help me insert and format a figure in LaTeX: $context"
            echo "Include proper placement, caption, label, and scaling options."
            ;;
        "bibliography")
            echo "Help me set up bibliography for: $context"
            echo "Include BibTeX setup, citation styles, and example entries."
            ;;
        "debug")
            echo "Debug this LaTeX error: $context"
            echo "Explain the error, why it occurs, and provide the fix."
            ;;
        *)
            echo "General LaTeX help for: $context"
            ;;
    esac
}

# 函数：LaTeX 公式助手
latex_formula() {
    echo -e "${CYAN}=== LaTeX 公式助手 ===${NC}"
    echo "常用公式类型："
    echo "1) 数学方程 (equation)"
    echo "2) 矩阵 (matrix)"
    echo "3) 积分/求和 (integral/sum)"
    echo "4) 分段函数 (piecewise)"
    echo "5) 定理/证明 (theorem/proof)"
    
    read -p "选择类型 (1-5): " choice
    read -p "描述你需要的公式: " description
    
    # 生成对应的 LaTeX 模板
    case $choice in
        1)
            cat << 'EOF'
% 数学方程示例
\begin{equation}
    f(x) = ax^2 + bx + c
    \label{eq:quadratic}
\end{equation}

% 多行对齐方程
\begin{align}
    f(x) &= (x+a)(x+b) \\
         &= x^2 + (a+b)x + ab
\end{align}
EOF
            ;;
        2)
            cat << 'EOF'
% 矩阵示例
\begin{equation}
    A = \begin{bmatrix}
        a_{11} & a_{12} & \cdots & a_{1n} \\
        a_{21} & a_{22} & \cdots & a_{2n} \\
        \vdots & \vdots & \ddots & \vdots \\
        a_{m1} & a_{m2} & \cdots & a_{mn}
    \end{bmatrix}
\end{equation}

% 需要包: \usepackage{amsmath}
EOF
            ;;
        3)
            cat << 'EOF'
% 积分示例
\begin{equation}
    \int_{a}^{b} f(x)\,dx = F(b) - F(a)
\end{equation}

% 求和示例
\begin{equation}
    \sum_{i=1}^{n} i = \frac{n(n+1)}{2}
\end{equation}
EOF
            ;;
        4)
            cat << 'EOF'
% 分段函数示例
\begin{equation}
    f(x) = \begin{cases}
        x^2 & \text{if } x \geq 0 \\
        -x^2 & \text{if } x < 0
    \end{cases}
\end{equation}

% 需要包: \usepackage{amsmath}
EOF
            ;;
        5)
            cat << 'EOF'
% 定理环境
\newtheorem{theorem}{Theorem}[section]
\newtheorem{lemma}[theorem]{Lemma}

\begin{theorem}[Fermat's Last Theorem]
    For $n > 2$, there are no three positive integers $a$, $b$, $c$ 
    that satisfy $a^n + b^n = c^n$.
\end{theorem}

\begin{proof}
    The proof is left as an exercise to the reader.
\end{proof}

% 需要包: \usepackage{amsthm}
EOF
            ;;
    esac
    
    echo -e "\n${GREEN}提示: 使用 Claude Code 询问具体的公式写法${NC}"
    generate_latex_prompt "formula" "$description"
}

# 函数：LaTeX 表格生成器
latex_table() {
    echo -e "${CYAN}=== LaTeX 表格生成器 ===${NC}"
    read -p "表格行数: " rows
    read -p "表格列数: " cols
    read -p "是否需要表头? (y/n): " header
    
    echo -e "\n${YELLOW}生成的表格模板:${NC}"
    
    # 生成列格式
    col_format="{"
    for ((i=1; i<=$cols; i++)); do
        col_format="${col_format}c"
        if [ $i -lt $cols ]; then
            col_format="${col_format}|"
        fi
    done
    col_format="${col_format}}"
    
    # 生成表格
    cat << EOF
\begin{table}[h]
    \centering
    \caption{Your Table Caption}
    \label{tab:example}
    \begin{tabular}${col_format}
        \hline
EOF
    
    if [ "$header" = "y" ]; then
        for ((i=1; i<=$cols; i++)); do
            printf "        Header $i"
            if [ $i -lt $cols ]; then
                printf " & "
            else
                printf " \\\\\\\\"
            fi
        done
        echo -e "\n        \\hline"
    fi
    
    for ((r=1; r<=$rows; r++)); do
        printf "        "
        for ((c=1; c<=$cols; c++)); do
            printf "Data R${r}C${c}"
            if [ $c -lt $cols ]; then
                printf " & "
            else
                printf " \\\\\\\\"
            fi
        done
        echo ""
    done
    
    cat << EOF
        \hline
    \end{tabular}
\end{table}
EOF
}

# 函数：LaTeX 错误诊断
latex_diagnose() {
    echo -e "${CYAN}=== LaTeX 错误诊断工具 ===${NC}"
    
    # 查找最近的日志文件
    LOG_FILE=$(ls -t *.log 2>/dev/null | head -n1)
    
    if [ -z "$LOG_FILE" ]; then
        echo -e "${RED}未找到日志文件${NC}"
        return 1
    fi
    
    echo -e "${YELLOW}分析日志文件: $LOG_FILE${NC}\n"
    
    # 提取错误和警告
    echo -e "${RED}错误:${NC}"
    grep -n "^!" "$LOG_FILE" | head -10 || echo "没有发现错误"
    
    echo -e "\n${YELLOW}警告:${NC}"
    grep -n "Warning:" "$LOG_FILE" | head -10 || echo "没有发现警告"
    
    echo -e "\n${BLUE}未定义引用:${NC}"
    grep -n "undefined" "$LOG_FILE" | head -10 || echo "没有未定义引用"
    
    echo -e "\n${MAGENTA}缺失文件:${NC}"
    grep -n "File.*not found" "$LOG_FILE" | head -10 || echo "没有缺失文件"
}

# 函数：LaTeX 包管理助手
latex_packages() {
    echo -e "${CYAN}=== LaTeX 常用包参考 ===${NC}"
    
    cat << 'EOF'
% === 基础包 ===
\usepackage[utf8]{inputenc}      % 输入编码
\usepackage[T1]{fontenc}         % 字体编码
\usepackage{babel}               % 多语言支持

% === 数学包 ===
\usepackage{amsmath}             % AMS 数学环境
\usepackage{amssymb}             % AMS 符号
\usepackage{amsthm}              % 定理环境
\usepackage{mathtools}           % 数学工具扩展

% === 图形包 ===
\usepackage{graphicx}            % 插入图片
\usepackage{tikz}                % 绘图
\usepackage{pgfplots}            % 数据图表
\usepackage{subfigure}           % 子图

% === 表格包 ===
\usepackage{booktabs}            % 专业表格
\usepackage{multirow}            % 跨行单元格
\usepackage{longtable}           % 跨页表格
\usepackage{array}               % 表格扩展

% === 引用包 ===
\usepackage{hyperref}            % 超链接
\usepackage{cite}                % 引用管理
\usepackage{natbib}              % 自然科学引用

% === 代码包 ===
\usepackage{listings}            % 代码高亮
\usepackage{algorithm}           % 算法环境
\usepackage{algpseudocode}       % 伪代码

% === 版式包 ===
\usepackage{geometry}            % 页面布局
\usepackage{fancyhdr}            % 页眉页脚
\usepackage{setspace}            % 行距控制
EOF
    
    echo -e "\n${GREEN}提示: 根据需要选择相应的包${NC}"
}

# 函数：快速 Overleaf 同步
quick_sync() {
    echo -e "${CYAN}=== 快速 Overleaf 同步 ===${NC}"
    
    # 检查是否在 Git 仓库
    if ! git rev-parse --git-dir > /dev/null 2>&1; then
        echo -e "${RED}错误: 不在 Git 仓库中${NC}"
        echo "请先使用 overleaf-setup.sh 克隆项目"
        return 1
    fi
    
    # 显示当前状态
    echo -e "${YELLOW}当前 Git 状态:${NC}"
    git status --short
    
    # 询问操作
    echo -e "\n选择操作:"
    echo "1) 拉取 Overleaf 更新"
    echo "2) 推送到 Overleaf"
    echo "3) 双向同步"
    
    read -p "选择 (1-3): " choice
    
    case $choice in
        1)
            git pull origin main || git pull origin master
            echo -e "${GREEN}✓ 已拉取最新更新${NC}"
            ;;
        2)
            git add -A
            read -p "提交信息: " msg
            git commit -m "$msg"
            git push origin main || git push origin master
            echo -e "${GREEN}✓ 已推送到 Overleaf${NC}"
            ;;
        3)
            git pull origin main || git pull origin master
            git add -A
            if ! git diff --staged --quiet; then
                read -p "提交信息: " msg
                git commit -m "$msg"
                git push origin main || git push origin master
            fi
            echo -e "${GREEN}✓ 双向同步完成${NC}"
            ;;
    esac
}

# 主菜单
main() {
    echo -e "${MAGENTA}╔══════════════════════════════════════╗${NC}"
    echo -e "${MAGENTA}║   Claude Code LaTeX Helper Suite    ║${NC}"
    echo -e "${MAGENTA}╚══════════════════════════════════════╝${NC}"
    echo ""
    echo "1) LaTeX 公式助手"
    echo "2) LaTeX 表格生成器"
    echo "3) LaTeX 错误诊断"
    echo "4) LaTeX 包管理参考"
    echo "5) 快速 Overleaf 同步"
    echo "6) 生成 Claude 提示词"
    echo "7) 退出"
    echo ""
    read -p "选择功能 (1-7): " choice
    
    case $choice in
        1) latex_formula ;;
        2) latex_table ;;
        3) latex_diagnose ;;
        4) latex_packages ;;
        5) quick_sync ;;
        6) 
            read -p "任务类型 (formula/table/figure/bibliography/debug): " task
            read -p "具体需求: " context
            generate_latex_prompt "$task" "$context"
            ;;
        7) echo -e "${GREEN}再见!${NC}"; exit 0 ;;
        *) echo -e "${RED}无效选择${NC}"; main ;;
    esac
    
    echo -e "\n${YELLOW}按回车继续...${NC}"
    read
    main
}

# 如果直接运行脚本，显示主菜单
if [ "${BASH_SOURCE[0]}" == "${0}" ]; then
    main
fi