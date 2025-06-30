"""
Display utilities for user interface.
"""

def print_welcome_message():
    """Print the PandaDock welcome message."""
    
    welcome_ascii = r"""
╔════════════════════════════════════════════════════════════════════════════════╗
║    ██████╗  █████╗ ███╗   ██╗██████╗  █████╗ ██████╗  ██████╗  ██████╗██╗  ██╗ ║
║    ██╔══██╗██╔══██╗████╗  ██║██╔══██╗██╔══██╗██╔══██╗██╔═══██╗██╔════╝██║ ██╔╝ ║
║    ██████╔╝███████║██╔██╗ ██║██║  ██║███████║██║  ██║██║   ██║██║     █████╔╝  ║
║    ██╔═══╝ ██╔══██║██║╚██╗██║██║  ██║██╔══██║██║  ██║██║   ██║██║     ██╔═██╗  ║
║    ██║     ██║  ██║██║ ╚████║██████╔╝██║  ██║██████╔╝╚██████╔╝╚██████╗██║  ██╗ ║
║    ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═══╝╚═════╝ ╚═╝  ╚═╝╚═════╝  ╚═════╝  ╚═════╝╚═╝  ╚═╝ ║
║                                                                                ║
║               🐼 Physics-based Molecular Docking Tool 🚀                       ║
║               Clean • Modular • CPU/GPU Ready                                  ║
║                                                                                ║
╚════════════════════════════════════════════════════════════════════════════════╝
    """
    
    print(welcome_ascii)
    print("🔬 Refactored Architecture - Easy to Debug & Extend")
    print("💻 Unified CPU/GPU Support - Works on Any System")
    print("🧬 Robust Pose Generation - Handles Failures Gracefully")
    print("")


def print_progress(current: int, total: int, operation: str = "Processing"):
    """Print progress bar."""
    
    if total == 0:
        return
    
    percent = (current / total) * 100
    filled = int(percent / 2)  # 50 character progress bar
    bar = "█" * filled + "░" * (50 - filled)
    
    print(f"\r{operation}: |{bar}| {percent:.1f}% ({current}/{total})", end="", flush=True)
    
    if current == total:
        print()  # New line when complete