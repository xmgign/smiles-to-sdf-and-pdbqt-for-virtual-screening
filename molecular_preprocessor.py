#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
分子预处理工具 - 精简版
将SMILES转换为SDF和PDBQT格式，用于虚拟筛选
"""

import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import logging
from pathlib import Path
from tqdm import tqdm
import multiprocessing as mp
from functools import partial
import time
from typing import List, Tuple, Optional, Dict, Any

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class MolecularPreprocessor:
    """分子预处理器 - 精简版"""
    
    def __init__(self, output_dir="output", max_workers=None):
        """
        初始化分子预处理器
        
        Args:
            output_dir: 输出目录
            max_workers: 最大工作进程数
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.max_workers = max_workers or mp.cpu_count()
        
        # 创建子目录
        self.sdf_dir = self.output_dir / "sdf_files"
        self.pdbqt_dir = self.output_dir / "pdbqt_files"
        self.sdf_dir.mkdir(exist_ok=True)
        self.pdbqt_dir.mkdir(exist_ok=True)
        
        # 统计信息
        self.stats = {
            'total_processed': 0,
            'sdf_successful': 0,
            'sdf_failed': 0,
            'pdbqt_successful': 0,
            'pdbqt_failed': 0,
            'force_field_stats': {'MMFF_success': 0, 'UFF_success': 0, 'no_optimization': 0}
        }
        
        # 检查Meeko
        self.meeko_available = self._check_meeko()
        
        logger.info(f"分子预处理器初始化完成")
        logger.info(f"输出目录: {self.output_dir}")
        logger.info(f"使用 {self.max_workers} 个进程")
        logger.info(f"Meeko可用性: {self.meeko_available}")
    
    def _check_meeko(self) -> bool:
        """检查Meeko是否可用"""
        try:
            import meeko
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            test_mol = Chem.MolFromSmiles("CCO")
            if test_mol:
                test_mol = Chem.AddHs(test_mol)
                AllChem.EmbedMolecule(test_mol)
                mk_prep = meeko.MoleculePreparation()
                mk_prep.prepare(test_mol)
                return True
            return False
        except:
            return False
    
    def detect_smiles_column(self, df: pd.DataFrame) -> Optional[str]:
        """智能检测SMILES列"""
        possible_names = [
            'smiles', 'SMILES', 'Smiles',
            'qsar_ready_smiles', 'QSAR_READY_SMILES',
            'canonical_smiles', 'Canonical_SMILES',
            'structure', 'Structure', 'molecule', 'Molecule'
        ]
        
        # 检查明确列名
        for col in df.columns:
            if col in possible_names:
                return col
        
        # 内容检测
        for col in df.columns:
            if df[col].dtype == 'object':
                sample_values = df[col].dropna().head(10)
                smiles_like_count = 0
                
                for value in sample_values:
                    if isinstance(value, str):
                        if any(char in value for char in ['C', 'N', 'O', 'S', 'P', '(', ')', '=', '#', '[', ']']):
                            smiles_like_count += 1
                
                if smiles_like_count / len(sample_values) > 0.5:
                    return col
        
        return None
    
    def generate_missing_columns(self, df: pd.DataFrame, smiles_col: str) -> pd.DataFrame:
        """生成缺少的列"""
        new_df = df.copy()
        
        required_columns = {
            'DTXSID': lambda i: f'COMPOUND_{i+1:06d}',
            'PREFERRED_NAME': lambda i: f'Compound_{i+1}',
            'CASRN': lambda i: '',
            'IUPAC_NAME': lambda i: '',
            'MOLECULAR_FORMULA': lambda i: '',
            'AVERAGE_MASS': lambda i: '',
            'QSAR_READY_SMILES': lambda i: df.iloc[i][smiles_col]
        }
        
        for col, generator in required_columns.items():
            if col not in new_df.columns:
                new_df[col] = [generator(i) for i in range(len(df))]
        
        return new_df
    
    def smiles_to_3d_mol(self, smiles: str, compound_id: str) -> Tuple[bool, Any, str]:
        """SMILES转3D分子"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return False, f"无法从SMILES创建分子: {smiles}", "invalid_smiles"
            
            mol = Chem.AddHs(mol)
            
            # 生成3D坐标
            success = False
            for attempt in range(3):
                try:
                    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, 
                                        useBasicKnowledge=True, randomSeed=42+attempt)
                    success = True
                    break
                except:
                    continue
            
            if not success:
                return False, "3D坐标生成失败", "embedding_failed"
            
            # 力场优化
            force_field_used = "none"
            try:
                result = AllChem.MMFFOptimizeMolecule(mol)
                if result == 0:
                    force_field_used = "MMFF"
            except:
                try:
                    result = AllChem.UFFOptimizeMolecule(mol)
                    if result == 0:
                        force_field_used = "UFF"
                except:
                    force_field_used = "none"
            
            mol.SetProp("DTXSID", str(compound_id))
            mol.SetProp("SMILES", smiles)
            mol.SetProp("FORCE_FIELD", force_field_used)
            
            return True, mol, force_field_used
            
        except Exception as e:
            return False, f"处理失败: {str(e)}", "error"
    
    def sdf_to_pdbqt(self, sdf_path: Path) -> Tuple[bool, str]:
        """SDF转PDBQT"""
        if not self.meeko_available:
            return False, "Meeko不可用"
        
        try:
            import meeko
            from rdkit import Chem
            
            supplier = Chem.SDMolSupplier(str(sdf_path))
            mol = next(supplier)
            
            if mol is None:
                return False, "无法读取SDF文件"
            
            mol = Chem.AddHs(mol)
            
            mk_prep = meeko.MoleculePreparation()
            setups = mk_prep.prepare(mol)
            
            if not setups:
                return False, "Meeko准备失败"
            
            from meeko import PDBQTWriterLegacy
            pdbqt_result = PDBQTWriterLegacy.write_string(setups[0])
            
            if isinstance(pdbqt_result, tuple):
                pdbqt_string = pdbqt_result[0]
            else:
                pdbqt_string = pdbqt_result
            
            pdbqt_path = self.pdbqt_dir / f"{sdf_path.stem}.pdbqt"
            with open(pdbqt_path, 'w') as f:
                f.write(pdbqt_string)
            
            return True, str(pdbqt_path)
            
        except Exception as e:
            return False, f"转换失败: {str(e)}"
    
    def process_compounds(self, df: pd.DataFrame, batch_size: int = 1000, 
                         convert_pdbqt: bool = True) -> int:
        """处理化合物"""
        # 检测SMILES列
        smiles_col = self.detect_smiles_column(df)
        if smiles_col is None:
            logger.error("无法识别SMILES列")
            return 0
        
        logger.info(f"检测到SMILES列: {smiles_col}")
        
        # 生成缺少的列
        df = self.generate_missing_columns(df, smiles_col)
        
        # 清理数据
        df = df.dropna(subset=['QSAR_READY_SMILES'])
        df = df[df['QSAR_READY_SMILES'] != '']
        
        logger.info(f"开始处理 {len(df)} 个化合物")
        
        # 分批处理
        total_successful = 0
        all_sdf_files = []
        
        for i in tqdm(range(0, len(df), batch_size), desc="处理批次"):
            batch_end = min(i + batch_size, len(df))
            batch_df = df.iloc[i:batch_end]
            
            # 处理批次
            with mp.Pool(processes=self.max_workers) as pool:
                results = list(tqdm(
                    pool.imap(self._process_single_compound, batch_df.iterrows()),
                    total=len(batch_df),
                    desc=f"批次 {i//batch_size + 1}"
                ))
            
            # 保存SDF文件并统计力场使用
            batch_successful = 0
            for result in results:
                if result['success']:
                    sdf_path = self.sdf_dir / f"{result['compound_id']}.sdf"
                    writer = Chem.SDWriter(str(sdf_path))
                    writer.write(result['mol'])
                    writer.close()
                    all_sdf_files.append(sdf_path)
                    batch_successful += 1
                    
                    # 统计力场使用情况
                    force_field = result['force_field_used']
                    if force_field == 'MMFF':
                        self.stats['force_field_stats']['MMFF_success'] += 1
                    elif force_field == 'UFF':
                        self.stats['force_field_stats']['UFF_success'] += 1
                    else:
                        self.stats['force_field_stats']['no_optimization'] += 1
            
            total_successful += batch_successful
            
            # 更新统计
            self.stats['total_processed'] += len(batch_df)
            self.stats['sdf_successful'] += batch_successful
            self.stats['sdf_failed'] += len(batch_df) - batch_successful
        
        logger.info(f"SDF生成完成！成功: {total_successful}")
        
        # 转换PDBQT
        if convert_pdbqt and all_sdf_files:
            logger.info("开始转换PDBQT文件...")
            pdbqt_successful = 0
            
            for sdf_path in tqdm(all_sdf_files, desc="SDF→PDBQT"):
                success, result = self.sdf_to_pdbqt(sdf_path)
                if success:
                    pdbqt_successful += 1
                else:
                    logger.debug(f"PDBQT转换失败 {sdf_path.name}: {result}")
            
            self.stats['pdbqt_successful'] = pdbqt_successful
            self.stats['pdbqt_failed'] = len(all_sdf_files) - pdbqt_successful
            logger.info(f"PDBQT转换完成！成功: {pdbqt_successful}")
        
        # 生成文件列表
        self._generate_file_lists()
        
        return total_successful
    
    def _process_single_compound(self, row_data):
        """处理单个化合物"""
        index, row = row_data
        compound_id = row['DTXSID']
        smiles = row['QSAR_READY_SMILES']
        
        success, mol_or_error, force_field_used = self.smiles_to_3d_mol(smiles, compound_id)
        
        return {
            'compound_id': compound_id,
            'success': success,
            'mol': mol_or_error if success else None,
            'error': mol_or_error if not success else None,
            'force_field_used': force_field_used
        }
    
    def _generate_file_lists(self):
        """生成文件列表"""
        # SDF文件列表
        sdf_files = list(self.sdf_dir.glob("*.sdf"))
        sdf_list_file = self.output_dir / "sdf_file_list.txt"
        with open(sdf_list_file, 'w') as f:
            for sdf_file in sdf_files:
                f.write(f"{sdf_file.name}\n")
        
        # PDBQT文件列表
        if self.meeko_available:
            pdbqt_files = list(self.pdbqt_dir.glob("*.pdbqt"))
            pdbqt_list_file = self.output_dir / "pdbqt_file_list.txt"
            with open(pdbqt_list_file, 'w') as f:
                for pdbqt_file in pdbqt_files:
                    f.write(f"{pdbqt_file.name}\n")
        
        logger.info(f"文件列表已生成: {len(sdf_files)} 个SDF文件")
    
    def generate_report(self):
        """生成处理报告"""
        report = f"""
分子预处理报告
==============

处理统计:
- 总处理化合物数: {self.stats['total_processed']:,}
- 成功生成SDF: {self.stats['sdf_successful']:,}
- SDF失败数量: {self.stats['sdf_failed']:,}
- 成功率: {self.stats['sdf_successful']/max(self.stats['total_processed'], 1)*100:.2f}%

力场优化统计:
- MMFF成功: {self.stats['force_field_stats']['MMFF_success']:,}
- UFF成功: {self.stats['force_field_stats']['UFF_success']:,}
- 无优化: {self.stats['force_field_stats']['no_optimization']:,}

PDBQT转换统计:
- 成功转换: {self.stats['pdbqt_successful']:,}
- 转换失败: {self.stats['pdbqt_failed']:,}
- Meeko可用性: {self.meeko_available}

输出文件:
- SDF文件: {self.sdf_dir}
- PDBQT文件: {self.pdbqt_dir}
- 文件列表: {self.output_dir}/sdf_file_list.txt
"""
        
        print(report)
        
        # 保存报告
        with open(self.output_dir / "processing_report.txt", "w", encoding="utf-8") as f:
            f.write(report)
        
        return report

def main():
    """主函数"""
    print("分子预处理工具 - 精简版")
    print("=" * 50)
    
    # 获取输入文件
    input_file = input("请输入SMILES文件路径: ").strip().strip('"')
    
    if not os.path.exists(input_file):
        print("文件不存在！")
        return
    
    # 加载数据
    try:
        if input_file.endswith('.csv'):
            df = pd.read_csv(input_file)
        elif input_file.endswith(('.xlsx', '.xls')):
            df = pd.read_excel(input_file)
        else:
            print("不支持的文件格式！")
            return
    except Exception as e:
        print(f"读取文件失败: {e}")
        return
    
    print(f"成功读取 {len(df)} 行数据")
    
    # 创建预处理器
    processor = MolecularPreprocessor()
    
    # 处理化合物
    start_time = time.time()
    successful_count = processor.process_compounds(df, convert_pdbqt=True)
    end_time = time.time()
    
    # 生成报告
    processor.generate_report()
    
    print(f"\n处理完成！")
    print(f"成功处理 {successful_count} 个化合物")
    print(f"耗时: {end_time - start_time:.2f} 秒")

if __name__ == "__main__":
    main()
