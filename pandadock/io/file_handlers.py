"""
File handling utilities for PandaDock.
"""

import logging
from pathlib import Path
from typing import Dict, Any, List


class FileHandlers:
    """Handles file I/O operations."""
    
    def __init__(self):
        """Initialize file handlers."""
        self.logger = logging.getLogger(__name__)
    
    def validate_input_files(self, protein_path: str, ligand_path: str) -> Dict[str, Any]:
        """
        Validate input files exist and are readable.
        
        Args:
            protein_path: Path to protein file
            ligand_path: Path to ligand file
            
        Returns:
            Validation result dictionary
        """
        
        result = {'valid': True, 'errors': []}
        
        # Check protein file
        if not Path(protein_path).exists():
            result['valid'] = False
            result['errors'].append(f"Protein file not found: {protein_path}")
        
        # Check ligand file  
        if not Path(ligand_path).exists():
            result['valid'] = False
            result['errors'].append(f"Ligand file not found: {ligand_path}")
        
        return result
    
    def create_output_directory(self, output_path: str) -> str:
        """
        Create output directory if it doesn't exist.
        
        Args:
            output_path: Path to output directory
            
        Returns:
            Absolute path to created directory
        """
        
        output_dir = Path(output_path)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        return str(output_dir.absolute())