# -*- coding: utf-8 -*-
"""
Modified on 2026-02-16
功能：处理 final 目录下的所有 *_HGT_statistics.txt 文件，
      过滤掉 recipient 和 donor 物种相同或为人类的记录，
      结果保存到 filter 目录。
"""

import os
import pandas as pd
from pathlib import Path

def extract_species_name(full_name):
    """
    从完整的物种描述中提取物种名称（属+种）
    例如: "Bifidobacterium pseudocatenulatum DSM 20438 = JCM 1200 = LMG 10505" -> "Bifidobacterium pseudocatenulatum"
    """
    if pd.isna(full_name) or full_name == "-" or full_name == "":
        return ""
    
    parts = str(full_name).strip().split()
    if len(parts) >= 2:
        # 检查第二个部分是否为"sp."或类似，这表示未鉴定到种
        if parts[1] in ["sp.", "sp", "spp.", "spp"]:
            return parts[0] + " " + parts[1] + ((" " + parts[2]) if len(parts) > 2 else "")
        else:
            return parts[0] + " " + parts[1]
    else:
        return str(full_name).strip()

def is_human_species(species_name):
    """
    检查物种名称是否包含人类(Homo sapiens)
    """
    if pd.isna(species_name) or species_name == "-" or species_name == "":
        return False
    
    species_lower = str(species_name).lower().strip()
    
    human_patterns = [
        "homo sapiens",
        "homo_sapiens",
        "h. sapiens",
        "h.sapiens",
        "human",
        "homo"
    ]
    
    for pattern in human_patterns:
        if pattern in species_lower:
            return True
    
    if species_lower.startswith("homo "):
        return True
    
    return False

def filter_hgt_files(input_base_dir, output_base_dir):
    """
    遍历输入目录，处理所有 *_HGT_statistics.txt 文件，过滤并保存到输出目录
    """
    input_base_path = Path(input_base_dir)
    output_base_path = Path(output_base_dir)
    output_base_path.mkdir(parents=True, exist_ok=True)
    
    total_files = 0
    processed_files = 0
    filtered_records = 0
    kept_records = 0
    human_records = 0
    
    # 递归遍历输入目录下的所有文件
    for root, dirs, files in os.walk(input_base_path):
        root_path = Path(root)
        
        for file_name in files:
            # 只处理以 _HGT_statistics.txt 结尾的文件
            if file_name.endswith('_HGT_statistics.txt'):
                input_file_path = root_path / file_name
                total_files += 1
                
                try:
                    # 读取文件
                    df = pd.read_csv(input_file_path, sep='\t')
                    
                    if df.empty:
                        print(f"文件为空: {input_file_path}")
                        continue
                    
                    # 检查必要的列是否存在
                    required_columns = ['recipient_species', 'donor_species']
                    if not all(col in df.columns for col in required_columns):
                        print(f"文件缺少必要列: {input_file_path}")
                        continue
                    
                    # 复制原始数据框用于处理
                    df_filtered = df.copy()
                    
                    # 提取物种名称（只到物种级别）
                    df_filtered['recipient_species_simple'] = df_filtered['recipient_species'].apply(extract_species_name)
                    df_filtered['donor_species_simple'] = df_filtered['donor_species'].apply(extract_species_name)
                    
                    # 创建过滤条件：排除完全相同的物种和物种级别相同的
                    condition_same_exact = df_filtered['recipient_species'] == df_filtered['donor_species']
                    condition_same_species = df_filtered['recipient_species_simple'] == df_filtered['donor_species_simple']
                    
                    # 检查是否为人类
                    condition_human_recipient = df_filtered['recipient_species'].apply(is_human_species)
                    condition_human_donor = df_filtered['donor_species'].apply(is_human_species)
                    condition_human = condition_human_recipient | condition_human_donor
                    
                    # 排除空字符串或只有"-"的情况
                    condition_valid_recipient = ~df_filtered['recipient_species_simple'].isin(["", "-", None])
                    condition_valid_donor = ~df_filtered['donor_species_simple'].isin(["", "-", None])
                    
                    # 应用过滤：保留不同物种且非人类的记录
                    df_filtered = df_filtered[
                        (~condition_same_exact) & 
                        (~condition_same_species) & 
                        (~condition_human) &
                        condition_valid_recipient & 
                        condition_valid_donor
                    ].copy()
                    
                    # 删除辅助列
                    if 'recipient_species_simple' in df_filtered.columns:
                        df_filtered = df_filtered.drop(['recipient_species_simple', 'donor_species_simple'], axis=1)
                    
                    # 统计信息
                    original_count = len(df)
                    filtered_count = len(df_filtered)
                    human_count = condition_human.sum()
                    
                    filtered_records += (original_count - filtered_count - human_count)
                    human_records += human_count
                    kept_records += filtered_count
                    
                    # 如果过滤后还有数据，保存到新文件
                    if not df_filtered.empty:
                        # 构建输出路径，保持原有目录结构
                        relative_path = input_file_path.relative_to(input_base_path)
                        output_file_path = output_base_path / relative_path
                        
                        # 确保输出目录存在
                        output_file_path.parent.mkdir(parents=True, exist_ok=True)
                        
                        # 保存过滤后的数据
                        df_filtered.to_csv(output_file_path, sep='\t', index=False)
                        processed_files += 1
                        print(f"已处理: {input_file_path} -> {output_file_path}")
                        print(f"  原始记录: {original_count}, 过滤后: {filtered_count}, 人类记录: {human_count}")
                    else:
                        print(f"文件过滤后为空，跳过: {input_file_path}")
                        if human_count > 0:
                            print(f"  其中包含人类记录: {human_count}")
                        
                except Exception as e:
                    print(f"处理文件时出错 {input_file_path}: {e}")
                    import traceback
                    traceback.print_exc()
    
    # 打印总结
    print("\n" + "="*50)
    print("处理完成!")
    print(f"总文件数: {total_files}")
    print(f"成功处理文件数: {processed_files}")
    print(f"总原始记录数: {filtered_records + kept_records + human_records}")
    print(f"过滤掉的记录数 (相同物种): {filtered_records}")
    print(f"过滤掉的记录数 (人类): {human_records}")
    print(f"保留的记录数: {kept_records}")
    print("="*50)

def main():
    # 设置路径
    input_base_dir = r"final"
    output_base_dir = r"filter"
    
    # 执行过滤
    filter_hgt_files(input_base_dir, output_base_dir)
    
    print(f"\n所有文件已处理完成!")
    print(f"过滤后的结果保存在: {output_base_dir}")

if __name__ == "__main__":
    main()