# -*- coding: utf-8 -*-
"""
Version management and update notification for PandaDock
"""

import requests
import sys
import json
import os
import time
from pathlib import Path
from typing import Optional, Dict, Any
import logging

__version__ = "2.5.0"

logger = logging.getLogger(__name__)


def get_version() -> str:
    """Get current PandaDock version"""
    return __version__


def check_for_updates(timeout: int = 5) -> Optional[Dict[str, Any]]:
    """
    Check for new PandaDock version on PyPI
    
    Args:
        timeout: Request timeout in seconds
        
    Returns:
        Dict with update info if available, None otherwise
    """
    try:
        response = requests.get(
            "https://pypi.org/pypi/pandadock/json",
            timeout=timeout,
            headers={'User-Agent': f'PandaDock/{__version__}'}
        )
        response.raise_for_status()
        
        data = response.json()
        latest_version = data['info']['version']
        
        if latest_version != __version__:
            return {
                'current_version': __version__,
                'latest_version': latest_version,
                'release_date': data['releases'][latest_version][0]['upload_time_iso_8601'],
                'summary': data['info']['summary'],
                'home_page': data['info']['home_page']
            }
    except Exception as e:
        logger.debug(f"Failed to check for updates: {e}")
    
    return None


def should_check_for_updates() -> bool:
    """
    Check if we should check for updates (max once per day)
    """
    try:
        update_file = Path.home() / '.pandadock' / 'last_update_check'
        if not update_file.exists():
            return True
            
        last_check = update_file.stat().st_mtime
        return (time.time() - last_check) > 86400  # 24 hours
    except Exception:
        return True


def save_update_check_time():
    """Save the last update check time"""
    try:
        update_dir = Path.home() / '.pandadock'
        update_dir.mkdir(exist_ok=True)
        update_file = update_dir / 'last_update_check'
        update_file.touch()
    except Exception as e:
        logger.debug(f"Failed to save update check time: {e}")


def get_update_message(update_info: Dict[str, Any]) -> str:
    """
    Generate user-friendly update notification message
    """
    current = update_info['current_version']
    latest = update_info['latest_version']
    
    message = f"""
╭─────────────────────────────────────────────────────────────╮
│                    🐼 PandaDock Update Available             │
├─────────────────────────────────────────────────────────────┤
│  Current version: {current:<20} Latest: {latest:<15} │
│                                                             │
│  🚀 New version available with improved performance         │
│     and better results! Update now for the best experience │
│                                                             │
│  📦 Update command:                                         │
│     pip install --upgrade pandadock                        │
│                                                             │
│  📖 Release info: {update_info.get('home_page', 'N/A'):<35} │
╰─────────────────────────────────────────────────────────────╯
"""
    return message


def check_and_notify_updates(silent: bool = False) -> bool:
    """
    Check for updates and notify user if available
    
    Args:
        silent: If True, don't print notification
        
    Returns:
        True if update is available, False otherwise
    """
    if not should_check_for_updates():
        return False
    
    update_info = check_for_updates()
    save_update_check_time()
    
    if update_info and not silent:
        print(get_update_message(update_info), file=sys.stderr)
        return True
    
    return update_info is not None


def print_version_info():
    """Print detailed version information"""
    print(f"PandaDock {__version__}")
    print(f"Python {sys.version}")
    print(f"Platform: {sys.platform}")
    
    # Check for updates
    update_info = check_for_updates()
    if update_info:
        print(f"\n⚠️  Update available: {update_info['latest_version']}")
        print("   Run 'pip install --upgrade pandadock' to update")
    else:
        print("\n✅ You're running the latest version!")