"""
Logging utilities for PandaDock.
"""

import logging
import os
from pathlib import Path
from datetime import datetime


def setup_logging(output_dir: str, log_level: str = "INFO") -> logging.Logger:
    """
    Set up logging for PandaDock.
    
    Args:
        output_dir: Directory for log files
        log_level: Logging level
        
    Returns:
        Configured logger instance
    """
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Create log file path
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = Path(output_dir) / f"pandadock_{timestamp}.log"
    
    # Configure logging
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()  # Also log to console
        ]
    )
    
    logger = logging.getLogger('pandadock')
    logger.info(f"Logging initialized - log file: {log_file}")
    
    return logger