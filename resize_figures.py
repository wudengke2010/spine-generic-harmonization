"""
一键修复所有图表 figsize，适配 SCI 期刊标准（全页宽 ≤7.0", 单栏 ≤3.35"）。

用法:
    cd E:/qm_harmonization_paper/derivatives/all-files
    python C:/Users/admin/WorkBuddy/2026-07-06-12-34-04/resize_figures.py

然后逐个运行修改后的脚本生成新图。
"""

import re
import os
import shutil
from pathlib import Path

TARGET_DIR = Path("E:/qm_harmonization_paper/derivatives/all-files")

# ---- 修改清单 ----
# (文件名, 行号, 旧figsize, 新figsize, 说明)
PATCHES = [
    # Fig1: 研究设计流程图 (竖版示意)
    # 当前 8.4×11.4" → 改为全页 6.5×8.8" (保持比例, 适合竖版)
    (
        "step5_8_fig1_schematic_v3.py", 303,
        "figsize=(8.4, 11.4)",
        "figsize=(6.5, 8.8)",
        "研究设计流程图 → 全页宽度 6.5\""
    ),

    # Fig2: 8个biomarker基线效应 (2行箱线+1行柱状, 最复杂)
    # 当前 10.0×8.0" → 改为全页 7.0×5.6"
    (
        "step5_7_fig2_baseline.py", 222,
        "figsize=(10.0, 8.0)",
        "figsize=(7.0, 5.6)",
        "基线效应多面板图 → 全页宽度 7.0\""
    ),

    # Fig3: 协调后性能 (2×2子图: log η² + ICC + 箱线 + PERMANOVA)
    # 当前 8.6×6.4" → 改为全页 7.0×5.2"
    (
        "step5_5_fig3_performance.py", 202,
        "figsize=(8.6, 6.4)",
        "figsize=(7.0, 5.2)",
        "性能评估多面板 → 全页宽度 7.0\""
    ),

    # Fig4: Pareto tradeoff 散点图 (2个cohort并列)
    # 当前 8.4×4.4" → 可改为单栏 3.35×1.8" 或全页 6.5×3.4"
    (
        "step5_4_fig4_tradeoff.py", 267,
        "figsize=(8.4, 4.4)",
        "figsize=(6.5, 3.4)",
        "Pareto前沿散点 → 全页宽度 6.5\""
    ),

    # Fig5: r_pre,post 热图 (2个cohort上下排列 + marginals)
    # 当前 9.0×8.6" → 改为全页 7.0×6.7"
    (
        "step5_6_fig5_heatmap_marginals.py", 231,
        "figsize=(9.0, 8.6)",
        "figsize=(7.0, 6.7)",
        "相关热图+边际分布 → 全页宽度 7.0\""
    ),

    # FigS: Design C稳健性 (1×3子图)
    # 当前 11.5×3.8" → 改为全页 7.0×2.3"
    (
        "step6_design_C_balanced.py", 374,
        "figsize=(11.5, 3.8)",
        "figsize=(7.0, 2.3)",
        "稳健性分析补充图 → 全页宽度 7.0\""
    ),
]


def patch_figsize(filepath: Path, line_num: int, old: str, new: str) -> bool:
    """替换指定行中的 figsize 值"""
    lines = filepath.read_text(encoding="utf-8").splitlines()
    if line_num < 1 or line_num > len(lines):
        print(f"  ✗ 行号 {line_num} 超出范围 (共 {len(lines)} 行)")
        return False

    target_line = lines[line_num - 1]
    if old not in target_line:
        print(f"  ✗ 第 {line_num} 行未找到 '{old}'")
        print(f"    实际内容: {target_line.strip()[:80]}")
        return False

    lines[line_num - 1] = target_line.replace(old, new)
    filepath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return True


def add_fontsize_boost(filepath: Path) -> bool:
    """
    如果脚本中设置了 rcParams 或 font.size，适当增大以补偿 figsize 缩小。
    仅在 scripts 没有显式设置 font.size 时添加。
    """
    content = filepath.read_text(encoding="utf-8")
    if "font.size" in content:
        # 已有设置，不覆盖
        return False

    # 在 import matplotlib 行之后插入字体设置
    # 查找 plt.rcParams 或 matplotlib 导入位置
    insert_after = None
    lines = content.splitlines()
    for i, line in enumerate(lines):
        if "import matplotlib" in line or "from matplotlib" in line:
            insert_after = i + 1
            break

    if insert_after is None:
        return False

    boost = "\n# --- 补偿 figsize 缩小后的字体可读性 ---\nimport matplotlib as mpl\nmpl.rcParams.update({'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8})\n"
    lines.insert(insert_after, boost)
    filepath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return True


def main():
    print("=" * 60)
    print("SCI 期刊图表尺寸修复工具")
    print("=" * 60)
    print(f"目标目录: {TARGET_DIR}")
    print(f"标准: 全页宽 ≤7.0\", 单栏 ≤3.35\"\n")

    # 1. 备份
    backup_dir = TARGET_DIR / "_figsize_backup"
    backup_dir.mkdir(exist_ok=True)
    print(f"[备份] → {backup_dir.name}/\n")

    # 2. 逐个打补丁
    patched_files = []
    for fname, line_num, old, new, desc in PATCHES:
        fp = TARGET_DIR / fname
        if not fp.exists():
            print(f"⏭ 跳过: {fname} (文件不存在)")
            continue

        # 备份
        shutil.copy2(fp, backup_dir / fname)
        print(f"📄 {fname}")
        print(f"   行 {line_num}: {old}")
        print(f"   →          {new}")
        print(f"   ({desc})")

        if patch_figsize(fp, line_num, old, new):
            print(f"   ✓ 已修改\n")
            patched_files.append(fp)
        else:
            print(f"   ✗ 修改失败\n")

    # 3. 字体适配 (可选)
    print("\n[字体适配] 以下脚本添加了基础字体设置以补偿尺寸缩小:")
    for fp in patched_files:
        if add_fontsize_boost(fp):
            print(f"  ✓ {fp.name}")

    # 4. 总结
    print("\n" + "=" * 60)
    print("修改完成!")
    print(f"备份位置: {backup_dir}")
    print(f"\n下一步 — 逐个运行以下脚本重新生成图表:")
    for fname, _, _, _, _ in PATCHES:
        fp = TARGET_DIR / fname
        if fp.exists():
            print(f"  python {fname}")
    print("\n生成后新图将覆盖原有 PNG/PDF。")
    print("=" * 60)

    # 5. 还原指南
    print(f"\n如需还原: 将 {backup_dir.name}/ 中的文件复制回上级目录即可。")


if __name__ == "__main__":
    main()
