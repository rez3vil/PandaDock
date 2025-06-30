"""
Input validation utilities.
"""

import os
from pathlib import Path
from typing import List, Dict, Any, Optional

def validate_input_files(protein_path: str, ligand_path: str, 
                        reference_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Validate input file paths and formats.
    
    Returns:
        Dictionary with validation results
    """
    
    validation_result = {
        'valid': True,
        'errors': [],
        'warnings': []
    }
    
    # Check protein file
    if not os.path.exists(protein_path):
        validation_result['valid'] = False
        validation_result['errors'].append(f"Protein file not found: {protein_path}")
    elif not protein_path.lower().endswith('.pdb'):
        validation_result['warnings'].append(f"Unusual protein file extension: {Path(protein_path).suffix}")
    
    # Check ligand file
    if not os.path.exists(ligand_path):
        validation_result['valid'] = False
        validation_result['errors'].append(f"Ligand file not found: {ligand_path}")
    elif not any(ligand_path.lower().endswith(ext) for ext in ['.sdf', '.mol', '.pdb']):
        validation_result['warnings'].append(f"Unusual ligand file extension: {Path(ligand_path).suffix}")
    
    # Check reference file if provided
    if reference_path:
        if not os.path.exists(reference_path):
            validation_result['valid'] = False
            validation_result['errors'].append(f"Reference file not found: {reference_path}")
    
    return validation_result