#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分子预处理工具使用示例
"""

from molecular_preprocessor import MolecularPreprocessor
import pandas as pd

def example_basic_usage():
    """基本使用示例"""
    print("=== 基本使用示例 ===")
    
    # 创建测试数据
    test_data = pd.DataFrame({
        'compound_id': ['COMP_001', 'COMP_002', 'COMP_003', 'COMP_004'],
        'smiles': ['CCO', 'CC(=O)O', 'c1ccccc1', 'CC(C)O'],
        'name': ['Ethanol', 'Acetic acid', 'Benzene', 'Isopropanol']
    })
    
    # 创建处理器
    processor = MolecularPreprocessor(
        output_dir="example_output",
        max_workers=2  # 使用2个进程
    )
    
    # 处理化合物
    success_count = processor.process_compounds(
        test_data, 
        batch_size=2,  # 小批次处理
        convert_pdbqt=True  # 生成PDBQT文件
    )
    
    print(f"成功处理 {success_count} 个化合物")
    
    # 生成详细报告
    processor.generate_report()

def example_large_file():
    """大文件处理示例"""
    print("=== 大文件处理示例 ===")
    
    # 模拟大文件数据
    import random
    simple_smiles = ['CCO', 'CC(=O)O', 'c1ccccc1', 'CC(C)O', 'CCN', 'CCC']
    
    large_data = pd.DataFrame({
        'id': [f'COMP_{i:06d}' for i in range(100)],
        'smiles': [random.choice(simple_smiles) for _ in range(100)]
    })
    
    # 创建处理器，使用更多进程
    processor = MolecularPreprocessor(
        output_dir="large_example_output",
        max_workers=4
    )
    
    # 批量处理
    success_count = processor.process_compounds(
        large_data,
        batch_size=20,  # 较大的批次
        convert_pdbqt=False  # 跳过PDBQT转换以提高速度
    )
    
    print(f"大文件处理：成功 {success_count}/100 个化合物")

if __name__ == "__main__":
    # 运行示例
    example_basic_usage()
    print("\n" + "="*50 + "\n")
    example_large_file()
