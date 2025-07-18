#!/usr/bin/env python3
"""
Test script for the Streamlit PandaDock application
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def test_imports():
    """Test that all required imports work"""
    try:
        import streamlit as st
        print("✅ Streamlit import successful")
        
        import pandas as pd
        print("✅ Pandas import successful")
        
        import numpy as np
        print("✅ NumPy import successful")
        
        import plotly.graph_objects as go
        print("✅ Plotly import successful")
        
        # Test PandaDock imports
        from pandadock.config import PandaDockConfig
        print("✅ PandaDock config import successful")
        
        from pandadock.docking.physics_engine import PhysicsEngine
        from pandadock.docking.ml_engine import MLEngine
        from pandadock.docking.ga_engine import GAEngine
        print("✅ PandaDock engine imports successful")
        
        from pandadock.reports.html_report import HTMLReportGenerator
        print("✅ PandaDock report generator import successful")
        
        return True
        
    except ImportError as e:
        print(f"❌ Import error: {e}")
        return False

def test_demo_functionality():
    """Test demo functionality"""
    try:
        # Import the demo function from our app
        from app import create_demo_results, DemoPose
        
        # Test demo pose creation
        demo_pose = DemoPose("test_pose", 8.5, -10.2, 5.5, 0.85)
        print(f"✅ Demo pose created: {demo_pose.pose_id}")
        
        # Test demo results creation
        demo_results = create_demo_results()
        print(f"✅ Demo results created: {len(demo_results)} poses")
        
        # Test demo pose methods
        affinity = demo_results[0].get_binding_affinity()
        ic50 = demo_results[0].get_ic50()
        print(f"✅ Demo pose methods work: affinity={affinity:.2f}, IC50={ic50:.2f}")
        
        return True
        
    except Exception as e:
        print(f"❌ Demo functionality error: {e}")
        return False

def test_configuration():
    """Test configuration creation"""
    try:
        from pandadock.config import PandaDockConfig
        
        config = PandaDockConfig()
        print("✅ PandaDock configuration created")
        
        # Test some basic configuration
        config.docking.mode = "balanced"
        config.docking.num_poses = 10
        config.scoring.scoring_function = "pandaml"
        
        print(f"✅ Configuration set: mode={config.docking.mode}, poses={config.docking.num_poses}")
        
        return True
        
    except Exception as e:
        print(f"❌ Configuration error: {e}")
        return False

def main():
    """Main test function"""
    print("🧪 Testing Streamlit PandaDock Application")
    print("=" * 50)
    
    # Test imports
    print("\n📦 Testing imports...")
    imports_ok = test_imports()
    
    # Test demo functionality  
    print("\n🎮 Testing demo functionality...")
    demo_ok = test_demo_functionality()
    
    # Test configuration
    print("\n⚙️ Testing configuration...")
    config_ok = test_configuration()
    
    # Summary
    print("\n" + "=" * 50)
    print("📊 Test Summary:")
    print(f"   Imports: {'✅ PASS' if imports_ok else '❌ FAIL'}")
    print(f"   Demo: {'✅ PASS' if demo_ok else '❌ FAIL'}")
    print(f"   Config: {'✅ PASS' if config_ok else '❌ FAIL'}")
    
    if imports_ok and demo_ok and config_ok:
        print("\n🎉 All tests passed! The Streamlit app should work correctly.")
        print("🚀 Run the app with: streamlit run app.py")
    else:
        print("\n❌ Some tests failed. Please check the errors above.")
        
    return imports_ok and demo_ok and config_ok

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)