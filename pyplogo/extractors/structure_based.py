import os
import warnings
import tempfile
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, DSSP
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio import Align
from Bio.Align import substitution_matrices
from dataclasses import dataclass, field
from typing import List, Tuple

@dataclass
class StructureData:
    """Container for structure extraction results"""
    sequence: str
    secondary_structure: str
    source_file: str
    method: str
    disulfide_bonds: List[Tuple[int, int]] = field(default_factory=list)  # 新增二硫键字段

class StructureExtractor:
    """Extract secondary structure and disulfide bonds from PDB/CIF files using DSSP"""
    
    def __init__(self):
        """Initialize structure extractor"""
        self.ss_mapping = {
            'H': 'H', 'G': 'G', 'I': 'I',  # Helices
            'E': 'E', 'B': 'E',            # Strands
            'T': 'T', 'S': 'S',            # Turns and bends
            ' ': 'C', '-': 'C'             # Coil
        }
        
        # 设置字典文件路径
        base_dir = os.path.dirname(os.path.abspath(__file__))
        dict_dir = os.path.join(base_dir, '..', 'data', 'dssp_dicts')
        dict_dir = os.path.abspath(dict_dir)
        os.environ['LIBCIFPP_DATA_DIR'] = dict_dir
        print(f"设置 DSSP 字典路径: {dict_dir}")
    
    def _run_dssp(self, model, file_path):
        """运行 DSSP 并处理可能的错误"""
        try:
            # 使用 mkdssp
            return DSSP(model, file_path, dssp='mkdssp')
        except Exception as e:
            # 尝试替代方法
            try:
                # 尝试使用标准 dssp
                return DSSP(model, file_path)
            except Exception as e2:
                # 尝试直接调用 mkdssp 可执行文件
                try:
                    from Bio.PDB.DSSP import dssp_dict_from_pdb_file
                    return dssp_dict_from_pdb_file(file_path, DSSP='mkdssp')
                except Exception as e3:
                    try:
                        # 尝试使用标准 dssp
                        return dssp_dict_from_pdb_file(file_path)
                    except Exception as e4:
                        raise RuntimeError(f"DSSP failed to produce an output: {str(e4)}") from e4
    
    def from_pdb(self, file_path: str, chain_id: str = 'A', input_sequence: str = None) -> StructureData:
        """从PDB文件提取二级结构和二硫键"""
        try:
            # 转换为绝对路径
            file_path = os.path.abspath(file_path)
            
            # 检查文件是否需要修复
            with open(file_path, 'r') as f:
                first_line = f.readline()
            
            # 如果缺少HEADER，创建修复版本
            if not first_line.startswith('HEADER'):
                print("检测到PDB文件缺少HEADER记录，正在修复...")
                file_path = self._fix_pdb_header(file_path)

            with warnings.catch_warnings():
                warnings.simplefilter("ignore", PDBConstructionWarning)
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure('protein', file_path)
                model = structure[0]
                
                # 运行DSSP
                dssp_result = self._run_dssp(model, file_path)
                
                # 提取结构数据
                result = self._extract_from_dssp(dssp_result, chain_id, file_path)
                
                # 提取二硫键信息
                result.disulfide_bonds = self._extract_disulfide_bonds(model, chain_id)
                
                # 序列匹配
                if input_sequence:
                    matched_ss = self.match_sequence(
                        result.sequence, 
                        result.secondary_structure, 
                        input_sequence
                    )
                    return StructureData(
                        sequence=input_sequence,
                        secondary_structure=matched_ss,
                        source_file=file_path,
                        method="dssp (matched)",
                        disulfide_bonds=result.disulfide_bonds
                    )
                
                return result
                
        except Exception as e:
            raise ValueError(f"Error processing PDB file: {e}")
        
    def _fix_pdb_header(self, original_file):
        """为PDB文件添加缺失的HEADER记录"""
        import tempfile
        
        # 读取原始内容
        with open(original_file, 'r') as f:
            content = f.read()
        
        # 添加HEADER记录
        header = "HEADER    PROTEIN                           01-JAN-70   FIXED_PDB\n"
        fixed_content = header + content
        
        # 创建临时文件
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
        temp_file.write(fixed_content)
        temp_file.close()
        
        return temp_file.name
    
    def from_cif(self, file_path: str, chain_id: str = 'A', input_sequence: str = None) -> StructureData:
        """从 CIF 文件提取二级结构和二硫键，可选自定义序列"""
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", PDBConstructionWarning)
                parser = MMCIFParser(QUIET=True)
                structure = parser.get_structure('protein', file_path)
                model = structure[0]  # 使用第一个模型
                
                # 运行 DSSP
                dssp = self._run_dssp(model, file_path)
                
                # 提取结构数据
                result = self._extract_from_dssp(dssp, chain_id, file_path)
                
                # 提取二硫键信息
                result.disulfide_bonds = self._extract_disulfide_bonds(model, chain_id)
                
                # 如果提供了自定义序列，进行匹配
                if input_sequence:
                    matched_ss = self.match_sequence(
                        result.sequence, 
                        result.secondary_structure, 
                        input_sequence
                    )
                    return StructureData(
                        sequence=input_sequence,
                        secondary_structure=matched_ss,
                        source_file=file_path,
                        method="dssp (matched)",
                        disulfide_bonds=result.disulfide_bonds  # 保留二硫键信息
                    )
                
                return result
                
        except Exception as e:
            raise ValueError(f"Error processing CIF file: {e}")
    
    def _extract_from_dssp(self, dssp, chain_id: str, file_path: str) -> StructureData:
        """从 DSSP 对象提取数据"""
        sequence = []
        secondary_structure = []
        
        # 处理不同的 DSSP 返回类型
        if isinstance(dssp, dict):
            # 来自 dssp_dict_from_pdb_file 的返回类型
            for key, value in dssp.items():
                if key[0] == chain_id:
                    aa = value[1]
                    ss = value[2]
                    
                    # 跳过非标准氨基酸
                    if aa not in 'ACDEFGHIKLMNPQRSTVWY':
                        continue
                        
                    sequence.append(aa)
                    secondary_structure.append(self.ss_mapping.get(ss, 'C'))
        else:
            # 来自 DSSP 类的返回类型
            for key in dssp.keys():
                if key[0] == chain_id:
                    aa = dssp[key][1]
                    ss = dssp[key][2]
                    
                    # 跳过非标准氨基酸
                    if aa not in 'ACDEFGHIKLMNPQRSTVWY':
                        continue
                        
                    sequence.append(aa)
                    secondary_structure.append(self.ss_mapping.get(ss, 'C'))
        
        if not sequence:
            raise ValueError(f"No data extracted for chain {chain_id}. "
                           f"Check if chain exists and contains standard amino acids.")
        
        return StructureData(
            sequence=''.join(sequence),
            secondary_structure=''.join(secondary_structure),
            source_file=file_path,
            method="dssp"
        )
    
    def _extract_disulfide_bonds(self, model, chain_id: str) -> List[Tuple[int, int]]:
        """从结构中提取二硫键信息"""
        disulfide_bonds = []
        
        # 获取指定链
        chain = None
        for c in model:
            if c.id == chain_id:
                chain = c
                break
        
        if not chain:
            return disulfide_bonds
        
        # 提取半胱氨酸残基
        cysteines = []
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname == 'CYS' and self._is_amino_acid(residue):
                cysteines.append(residue)
        
        # 检测二硫键 (距离小于3.0Å)
        for i, cys1 in enumerate(cysteines):
            for j, cys2 in enumerate(cysteines):
                if i >= j:
                    continue
                
                # 计算SG原子间距离
                sg1 = None
                sg2 = None
                for atom in cys1:
                    if atom.get_name() == 'SG':
                        sg1 = atom
                for atom in cys2:
                    if atom.get_name() == 'SG':
                        sg2 = atom
                
                if sg1 and sg2:
                    distance = sg1 - sg2
                    if distance < 3.0:  # 二硫键典型距离
                        # 获取残基索引 (基于PDB文件中的序号)
                        res_id1 = cys1.get_id()[1]
                        res_id2 = cys2.get_id()[1]
                        disulfide_bonds.append((res_id1, res_id2))
        
        return disulfide_bonds
    
    def _is_amino_acid(self, residue):
        """检查残基是否是标准氨基酸"""
        try:
            return residue.get_resname() in [
                'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                'THR', 'TRP', 'TYR', 'VAL'
            ]
        except:
            return False
    
    def match_sequence(self, structure_sequence, structure_ss, input_sequence):
        """
        匹配输入序列与结构序列，生成对应的二级结构字符串
        
        Args:
            structure_sequence: 从结构文件中提取的序列
            structure_ss: 从结构文件中提取的二级结构
            input_sequence: 用户输入的序列
            
        Returns:
            与输入序列长度匹配的二级结构字符串
        """
        # 如果输入序列为空，直接返回结构序列
        if not input_sequence:
            return structure_ss
        
        # 使用局部序列比对找到最佳匹配
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
        alignments = aligner.align(structure_sequence, input_sequence)
        
        if not alignments:
            # 如果没有找到匹配，返回全coil结构
            return 'C' * len(input_sequence)
        
        # 取最佳比对结果
        best_alignment = alignments[0]
        
        # 提取比对信息
        aligned_target = best_alignment.target
        aligned_query = best_alignment.query
        
        # 创建与输入序列等长的二级结构数组，初始化为coil
        matched_ss = ['C'] * len(input_sequence)
        
        # 遍历比对结果，填充匹配位置的二级结构
        target_idx = 0
        query_idx = 0
        
        for t_char, q_char in zip(aligned_target, aligned_query):
            if t_char != '-' and q_char != '-':
                # 匹配位置，使用结构中的二级结构
                if target_idx < len(structure_ss):
                    matched_ss[query_idx] = structure_ss[target_idx]
                target_idx += 1
                query_idx += 1
            elif t_char == '-' and q_char != '-':
                # 插入位置，保持为coil
                query_idx += 1
            elif t_char != '-' and q_char == '-':
                # 删除位置，跳过
                target_idx += 1
        
        return ''.join(matched_ss)