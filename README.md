# 分子预处理工具

一个简单高效的分子预处理工具，用于将SMILES字符串转换为3D SDF和PDBQT文件。

## 🚀 快速开始

### 1. 下载文件
只需要下载 `molecular_preprocessor.py` 这一个文件即可！

### 2. 安装依赖
```bash
pip install rdkit pandas tqdm
# 可选：安装Meeko用于PDBQT转换
pip install meeko
```

### 3. 直接使用
```python
python molecular_preprocessor.py
```

程序会提示您输入SMILES文件路径，支持格式：
- `.xlsx` - Excel文件
- `.csv` - CSV文件  
- `.smi` - SMILES文件
- `.txt` - 文本文件

## 📁 输出文件

程序会在 `output/` 目录下生成：
- `sdf_files/` - 3D SDF结构文件
- `pdbqt_files/` - PDBQT格式文件（如果安装了Meeko）
- `sdf_file_list.txt` - SDF文件列表
- `pdbqt_file_list.txt` - PDBQT文件列表
- `processing_report.txt` - 处理报告

## ✨ 主要特性

- 🧠 **智能SMILES检测**：自动识别SMILES列
- 🔧 **力场优化**：MMFF → UFF 自动回退
- ⚡ **多进程处理**：充分利用CPU性能
- 📊 **详细统计**：力场使用、成功率等
- 🎯 **PDBQT转换**：支持分子对接格式
- 📝 **进度显示**：实时显示处理进度

## 🔬 使用示例

### 基本用法
```python
from molecular_preprocessor import MolecularPreprocessor
import pandas as pd

# 创建处理器
processor = MolecularPreprocessor(output_dir="my_output")

# 准备数据（包含SMILES列）
data = pd.DataFrame({
    'compound_id': ['1', '2', '3'],
    'smiles': ['CCO', 'CC(=O)O', 'c1ccccc1']
})

# 处理化合物
success_count = processor.process_compounds(data, convert_pdbqt=True)
print(f"成功处理 {success_count} 个化合物")

# 生成报告
processor.generate_report()
```

### 批量处理大文件
```python
# 处理大文件时使用较小的批次大小
processor = MolecularPreprocessor(max_workers=4)
success_count = processor.process_compounds(data, batch_size=100)
```

## 📋 系统要求

- Python 3.7+
- RDKit
- Pandas
- tqdm
- Meeko (可选，用于PDBQT转换)

## 🐛 常见问题

**Q: 为什么有些化合物处理失败？**
A: 可能原因：无效的SMILES、3D坐标生成失败、力场优化失败

**Q: 可以处理多大的文件？**
A: 理论上没有限制，程序使用批处理和进度保存，支持大文件处理

**Q: PDBQT转换失败怎么办？**
A: 确保安装了Meeko：`pip install meeko`

## 📈 性能参考

- 测试环境：8核CPU
- 处理速度：约60个化合物/秒
- 内存使用：约500MB（1000个化合物批次）

## 📄 许可证

MIT License - 自由使用和修改

---

**简单、高效、易用！只需一个Python文件即可开始分子预处理！** 🧪✨
