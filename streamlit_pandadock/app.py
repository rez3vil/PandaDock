#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Streamlit PandaDock: Web interface for PandaDock molecular docking
"""

import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import os
import sys
import json
import tempfile
import zipfile
from pathlib import Path
import logging
from typing import Dict, List, Optional, Any
import io
import base64

# Add parent directory to path for PandaDock imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from pandadock.config import PandaDockConfig
from pandadock.docking.physics_engine import PhysicsEngine
from pandadock.docking.ml_engine import MLEngine
from pandadock.docking.ga_engine import GAEngine
from pandadock.reports.html_report import HTMLReportGenerator

# Configure page
st.set_page_config(
    page_title="PandaDock - Molecular Docking Platform",
    page_icon="🐼",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
.main-header {
    font-size: 2.5rem;
    font-weight: bold;
    color: #2E86AB;
    text-align: center;
    margin-bottom: 2rem;
}
.sub-header {
    font-size: 1.5rem;
    font-weight: bold;
    color: #A23B72;
    margin-bottom: 1rem;
}
.metric-card {
    background-color: #f0f2f6;
    padding: 1rem;
    border-radius: 0.5rem;
    margin-bottom: 1rem;
}
.success-box {
    background-color: #d4edda;
    border: 1px solid #c3e6cb;
    color: #155724;
    padding: 1rem;
    border-radius: 0.25rem;
    margin: 1rem 0;
}
.error-box {
    background-color: #f8d7da;
    border: 1px solid #f5c6cb;
    color: #721c24;
    padding: 1rem;
    border-radius: 0.25rem;
    margin: 1rem 0;
}
</style>
""", unsafe_allow_html=True)

def initialize_session_state():
    """Initialize session state variables"""
    if 'docking_results' not in st.session_state:
        st.session_state.docking_results = None
    if 'config' not in st.session_state:
        st.session_state.config = None
    if 'engine' not in st.session_state:
        st.session_state.engine = None
    if 'generated_files' not in st.session_state:
        st.session_state.generated_files = {}
    if 'output_dir' not in st.session_state:
        st.session_state.output_dir = None

def create_header():
    """Create main header with logo and title"""
    st.markdown('<div class="main-header">🐼 PandaDock - Molecular Docking Platform</div>', unsafe_allow_html=True)
    st.markdown("**Multi-Strategy, High-Performance Molecular Docking with Advanced Visualization**")
    st.markdown("---")

def create_sidebar():
    """Create sidebar with configuration options"""
    st.sidebar.markdown("## 🔧 Configuration")
    
    # Docking Mode Selection
    st.sidebar.markdown("### Docking Algorithm")
    mode = st.sidebar.selectbox(
        "Select Docking Mode",
        ["balanced", "precise", "fast"],
        index=0,
        help="Balanced: PandaML (ML-based), Precise: PandaPhysics (physics-based), Fast: PandaCore (evolutionary)"
    )
    
    # Scoring Function
    scoring_mapping = {
        "balanced": "pandaml",
        "precise": "pandaphysics", 
        "fast": "pandacore"
    }
    scoring = st.sidebar.selectbox(
        "Scoring Function",
        ["pandaml", "pandaphysics", "pandacore"],
        index=["pandaml", "pandaphysics", "pandacore"].index(scoring_mapping[mode]),
        help="PandaML: ML-based scoring, PandaPhysics: Physics-based, PandaCore: Evolutionary"
    )
    
    # Docking Parameters
    st.sidebar.markdown("### Docking Parameters")
    num_poses = st.sidebar.slider("Number of Poses", 1, 50, 10, help="Number of poses to generate")
    exhaustiveness = st.sidebar.slider("Exhaustiveness", 1, 32, 8, help="Search exhaustiveness (higher = more thorough)")
    energy_range = st.sidebar.slider("Energy Range (kcal/mol)", 1.0, 10.0, 3.0, help="Energy range for pose selection")
    
    # Grid Box Parameters
    st.sidebar.markdown("### Grid Box Settings")
    use_auto_center = st.sidebar.checkbox("Auto-detect binding site", value=True)
    
    if not use_auto_center:
        col1, col2, col3 = st.sidebar.columns(3)
        with col1:
            center_x = st.number_input("Center X", value=0.0)
        with col2:
            center_y = st.number_input("Center Y", value=0.0)
        with col3:
            center_z = st.number_input("Center Z", value=0.0)
    else:
        center_x = center_y = center_z = None
    
    size_x = st.sidebar.slider("Box Size X (Å)", 10.0, 50.0, 22.5)
    size_y = st.sidebar.slider("Box Size Y (Å)", 10.0, 50.0, 22.5)
    size_z = st.sidebar.slider("Box Size Z (Å)", 10.0, 50.0, 22.5)
    
    # Flexible Docking
    st.sidebar.markdown("### Flexible Docking")
    flexible_residues = st.sidebar.text_input(
        "Flexible Residues",
        placeholder="e.g., HIS57,SER195,TYR191",
        help="Comma-separated list of flexible residues"
    )
    side_chain_flexibility = st.sidebar.checkbox("Enable side chain flexibility", value=False)
    
    # Advanced Options
    st.sidebar.markdown("### Advanced Options")
    ml_rescoring = st.sidebar.checkbox("ML Rescoring", value=False, help="Use ML-based rescoring")
    seed = st.sidebar.number_input("Random Seed", value=42, help="Random seed for reproducibility")
    
    return {
        'mode': mode,
        'scoring': scoring,
        'num_poses': num_poses,
        'exhaustiveness': exhaustiveness,
        'energy_range': energy_range,
        'center': (center_x, center_y, center_z) if not use_auto_center else None,
        'size': (size_x, size_y, size_z),
        'flexible_residues': flexible_residues,
        'side_chain_flexibility': side_chain_flexibility,
        'ml_rescoring': ml_rescoring,
        'seed': seed
    }

def create_file_upload_section():
    """Create file upload section"""
    st.markdown('<div class="sub-header">📁 Input Files</div>', unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**Protein Receptor (PDB)**")
        protein_file = st.file_uploader(
            "Upload protein receptor",
            type=['pdb'],
            help="PDB file of the protein receptor"
        )
        
        if protein_file:
            st.success(f"✅ Protein loaded: {protein_file.name}")
            
    with col2:
        st.markdown("**Ligand Input**")
        ligand_option = st.radio(
            "Input type:",
            ["Single ligand file", "SMILES string", "Ligand library"]
        )
        
        if ligand_option == "Single ligand file":
            ligand_file = st.file_uploader(
                "Upload ligand",
                type=['sdf', 'mol2', 'pdb'],
                help="SDF, MOL2, or PDB file of the ligand"
            )
            ligand_smiles = None
            ligand_library = None
            
        elif ligand_option == "SMILES string":
            ligand_smiles = st.text_input(
                "Enter SMILES:",
                placeholder="e.g., CCO for ethanol",
                help="SMILES string of the ligand"
            )
            ligand_file = None
            ligand_library = None
            
        else:  # Ligand library
            ligand_library = st.file_uploader(
                "Upload ligand library",
                type=['smi', 'txt', 'csv'],
                help="File containing multiple ligands (SMILES format)"
            )
            ligand_file = None
            ligand_smiles = None
    
    return protein_file, ligand_file, ligand_smiles, ligand_library

def create_config_from_params(params: Dict, protein_file, ligand_file, ligand_smiles, output_dir: str) -> PandaDockConfig:
    """Create PandaDock configuration from parameters"""
    config = PandaDockConfig()
    
    # Docking parameters
    config.docking.mode = params['mode']
    config.docking.num_poses = params['num_poses']
    config.docking.exhaustiveness = params['exhaustiveness']
    config.docking.energy_range = params['energy_range']
    config.docking.seed = params['seed']
    
    # Flexible docking
    if params['flexible_residues']:
        config.docking.flexible_residues = [res.strip() for res in params['flexible_residues'].split(',')]
    config.docking.side_chain_flexibility = params['side_chain_flexibility']
    
    # Scoring
    config.scoring.scoring_function = params['scoring']
    config.scoring.use_ml_rescoring = params['ml_rescoring']
    
    # Grid box
    if params['center']:
        config.io.center_x, config.io.center_y, config.io.center_z = params['center']
    config.io.size_x, config.io.size_y, config.io.size_z = params['size']
    
    # IO settings
    config.io.output_dir = output_dir
    config.io.save_poses = True
    config.io.save_complex = True
    config.io.report_format = 'html'
    
    # Enable all visualization options
    config.enable_plots = True
    config.enable_interaction_maps = True
    config.enable_master_plot = True
    config.enable_txt_report = True
    config.enable_pandamap = True
    config.enable_pandamap_3d = True
    config.comprehensive_output = True
    
    return config

def get_docking_engine(config: PandaDockConfig):
    """Get appropriate docking engine based on configuration"""
    if config.docking.mode == 'precise':
        return PhysicsEngine(config)
    elif config.docking.mode == 'balanced':
        return MLEngine(config)
    elif config.docking.mode == 'fast':
        return GAEngine(config)
    else:
        raise ValueError(f"Unknown docking mode: {config.docking.mode}")

class DemoPose:
    """Demo pose class for demonstration purposes"""
    def __init__(self, pose_id: str, score: float, binding_affinity: float, ic50: float, confidence: float):
        self.pose_id = pose_id
        self.score = score
        self.binding_affinity = binding_affinity
        self.ic50 = ic50
        self.confidence = confidence
    
    def get_binding_affinity(self):
        return self.binding_affinity
    
    def get_ic50(self, units='uM'):
        if units == 'uM':
            return self.ic50
        elif units == 'nM':
            return self.ic50 * 1000
        return self.ic50

def create_demo_results():
    """Create demo results for testing the interface"""
    demo_poses = []
    
    # Generate realistic demo data
    np.random.seed(42)
    for i in range(10):
        pose_id = f"pose_{i+1}"
        score = np.random.uniform(6.0, 9.5)
        binding_affinity = np.random.uniform(-12.0, -6.0)
        ic50 = np.random.uniform(0.01, 100.0)
        confidence = np.random.uniform(0.5, 0.95)
        
        demo_poses.append(DemoPose(pose_id, score, binding_affinity, ic50, confidence))
    
    # Sort by score (highest first)
    demo_poses.sort(key=lambda x: x.score, reverse=True)
    
    return demo_poses

def run_docking_analysis(config: PandaDockConfig, protein_path: str, ligand_path: str):
    """Run docking analysis"""
    try:
        # Initialize engine
        engine = get_docking_engine(config)
        
        # Run docking
        results = engine.dock(protein_path, ligand_path)
        
        # Save poses
        poses_dir = os.path.join(config.io.output_dir, "poses")
        os.makedirs(poses_dir, exist_ok=True)
        engine.save_poses(results, poses_dir)
        
        # Generate comprehensive report
        report_generator = HTMLReportGenerator(config)
        
        # Algorithm and command info
        algorithm_info = {
            'algorithm': 'PandaDock',
            'version': '1.0.0',
            'scoring_function': config.scoring.scoring_function.upper(),
            'engine': config.docking.mode.title(),
            'mode': config.docking.mode
        }
        
        command_info = {
            'command': 'Streamlit Web Interface',
            'protein': protein_path,
            'ligand': ligand_path,
            'center': f"{getattr(config.io, 'center_x', 'auto')} {getattr(config.io, 'center_y', 'auto')} {getattr(config.io, 'center_z', 'auto')}",
            'size': f"{config.io.size_x} {config.io.size_y} {config.io.size_z}",
            'exhaustiveness': config.docking.exhaustiveness
        }
        
        # Generate comprehensive report
        generated_files = report_generator.generate_comprehensive_report(
            results=results,
            output_dir=config.io.output_dir,
            protein_name=Path(protein_path).stem,
            ligand_name=Path(ligand_path).stem,
            algorithm_info=algorithm_info,
            command_info=command_info
        )
        
        return results, generated_files
        
    except Exception as e:
        st.error(f"Docking failed: {str(e)}")
        return None, None

def create_results_visualization(results: List, generated_files: Dict):
    """Create results visualization section"""
    if not results:
        return
    
    st.markdown('<div class="sub-header">📊 Results Analysis</div>', unsafe_allow_html=True)
    
    # Create results DataFrame
    results_data = []
    for i, pose in enumerate(results):
        binding_affinity = pose.get_binding_affinity()
        ic50_um = pose.get_ic50(units='uM')
        
        results_data.append({
            'Rank': i + 1,
            'Pose_ID': pose.pose_id,
            'Score': pose.score,
            'Binding_Affinity': binding_affinity,
            'IC50_uM': ic50_um if ic50_um != float('inf') else None,
            'Confidence': getattr(pose, 'confidence', 0.0)
        })
    
    df = pd.DataFrame(results_data)
    
    # Summary metrics
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Total Poses", len(results))
    
    with col2:
        best_score = df['Score'].max()
        st.metric("Best Score", f"{best_score:.3f}")
    
    with col3:
        best_affinity = df['Binding_Affinity'].min()
        st.metric("Best Affinity", f"{best_affinity:.2f} kcal/mol")
    
    with col4:
        avg_confidence = df['Confidence'].mean()
        st.metric("Avg Confidence", f"{avg_confidence:.3f}")
    
    # Interactive plots
    tab1, tab2, tab3, tab4 = st.tabs(["📈 Score Distribution", "🎯 Binding Affinity", "💊 IC50 Analysis", "📋 Results Table"])
    
    with tab1:
        # Score distribution plot
        fig_score = px.histogram(df, x='Score', nbins=20, title='Docking Score Distribution')
        fig_score.update_layout(xaxis_title='Docking Score', yaxis_title='Frequency')
        st.plotly_chart(fig_score, use_container_width=True)
        
        # Score vs Rank scatter plot
        fig_rank = px.scatter(df, x='Rank', y='Score', color='Confidence', 
                             title='Docking Score vs Pose Rank',
                             hover_data=['Pose_ID', 'Binding_Affinity'])
        st.plotly_chart(fig_rank, use_container_width=True)
    
    with tab2:
        # Binding affinity plot
        fig_affinity = px.bar(df.head(10), x='Pose_ID', y='Binding_Affinity', 
                             title='Top 10 Poses - Binding Affinity',
                             color='Binding_Affinity', color_continuous_scale='RdYlBu_r')
        fig_affinity.update_layout(xaxis_title='Pose ID', yaxis_title='Binding Affinity (kcal/mol)')
        st.plotly_chart(fig_affinity, use_container_width=True)
        
        # Affinity vs Score correlation
        fig_corr = px.scatter(df, x='Score', y='Binding_Affinity', 
                             title='Binding Affinity vs Docking Score',
                             hover_data=['Pose_ID', 'Rank'])
        st.plotly_chart(fig_corr, use_container_width=True)
    
    with tab3:
        # IC50 analysis (only for poses with finite IC50)
        df_ic50 = df[df['IC50_uM'].notna() & (df['IC50_uM'] != float('inf'))]
        
        if not df_ic50.empty:
            # IC50 distribution
            fig_ic50 = px.histogram(df_ic50, x='IC50_uM', nbins=20, title='IC50 Distribution')
            fig_ic50.update_layout(xaxis_title='IC50 (µM)', yaxis_title='Frequency')
            st.plotly_chart(fig_ic50, use_container_width=True)
            
            # IC50 vs Binding Affinity
            fig_ic50_affinity = px.scatter(df_ic50, x='Binding_Affinity', y='IC50_uM', 
                                          title='IC50 vs Binding Affinity',
                                          hover_data=['Pose_ID', 'Rank'])
            fig_ic50_affinity.update_layout(xaxis_title='Binding Affinity (kcal/mol)', yaxis_title='IC50 (µM)')
            st.plotly_chart(fig_ic50_affinity, use_container_width=True)
        else:
            st.info("No valid IC50 data available for visualization.")
    
    with tab4:
        # Results table
        st.dataframe(df, use_container_width=True)
        
        # Export options
        csv = df.to_csv(index=False)
        st.download_button(
            label="📥 Download Results as CSV",
            data=csv,
            file_name="pandadock_results.csv",
            mime="text/csv"
        )

def create_file_downloads(generated_files: Dict, output_dir: str):
    """Create file download section"""
    if not generated_files:
        return
    
    st.markdown('<div class="sub-header">📥 Download Results</div>', unsafe_allow_html=True)
    
    # Create a zip file with all results
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
        for file_type, file_path in generated_files.items():
            # Handle both single files and lists of files
            if isinstance(file_path, list):
                for fp in file_path:
                    if os.path.exists(fp):
                        zip_file.write(fp, os.path.basename(fp))
            else:
                if os.path.exists(file_path):
                    zip_file.write(file_path, os.path.basename(file_path))
        
        # Add poses directory if it exists
        poses_dir = os.path.join(output_dir, "poses")
        if os.path.exists(poses_dir):
            for root, dirs, files in os.walk(poses_dir):
                for file in files:
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, output_dir)
                    zip_file.write(file_path, arcname)
    
    zip_buffer.seek(0)
    
    st.download_button(
        label="📦 Download All Results (ZIP)",
        data=zip_buffer.getvalue(),
        file_name="pandadock_results.zip",
        mime="application/zip"
    )
    
    # Individual file downloads
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**📊 Analysis Files**")
        for file_type, file_path in generated_files.items():
            if 'plot' in file_type or 'analysis' in file_type:
                # Handle both single files and lists of files
                file_paths = file_path if isinstance(file_path, list) else [file_path]
                for fp in file_paths:
                    if os.path.exists(fp):
                        with open(fp, 'rb') as f:
                            file_data = f.read()
                        
                        file_ext = Path(fp).suffix
                        mime_type = 'image/png' if file_ext == '.png' else 'application/pdf' if file_ext == '.pdf' else 'text/html'
                        
                        st.download_button(
                            label=f"📥 {Path(fp).name}",
                            data=file_data,
                            file_name=Path(fp).name,
                            mime=mime_type
                        )
    
    with col2:
        st.markdown("**📄 Report Files**")
        for file_type, file_path in generated_files.items():
            if 'report' in file_type or file_type in ['html', 'json', 'csv']:
                # Handle both single files and lists of files
                file_paths = file_path if isinstance(file_path, list) else [file_path]
                for fp in file_paths:
                    if os.path.exists(fp):
                        with open(fp, 'rb') as f:
                            file_data = f.read()
                        
                        file_ext = Path(fp).suffix
                        mime_type = 'text/html' if file_ext == '.html' else 'application/json' if file_ext == '.json' else 'text/csv'
                        
                        st.download_button(
                            label=f"📥 {Path(fp).name}",
                            data=file_data,
                            file_name=Path(fp).name,
                            mime=mime_type
                        )

def display_generated_plots(generated_files: Dict):
    """Display generated plots in the interface"""
    if not generated_files:
        return
    
    st.markdown('<div class="sub-header">🎨 Generated Visualizations</div>', unsafe_allow_html=True)
    
    # Display master publication plot
    if 'master_publication' in generated_files:
        master_plot_path = generated_files['master_publication']
        # Handle both single files and lists of files
        if isinstance(master_plot_path, list):
            master_plot_path = master_plot_path[0] if master_plot_path else None
        
        if master_plot_path and os.path.exists(master_plot_path):
            st.markdown("**📊 Master Publication Plot**")
            st.image(master_plot_path, caption="Comprehensive Analysis Dashboard", use_container_width=True)
    
    # Display other plots in tabs
    plot_files = {}
    for k, v in generated_files.items():
        if 'plot' in k.lower() or 'analysis' in k.lower():
            # Handle both single files and lists of files
            if isinstance(v, list):
                plot_files[k] = v[0] if v else None
            else:
                plot_files[k] = v
    
    if plot_files:
        plot_tabs = st.tabs([f"📈 {k.replace('_', ' ').title()}" for k in plot_files.keys()])
        
        for i, (plot_type, plot_path) in enumerate(plot_files.items()):
            if plot_path and os.path.exists(plot_path) and plot_type != 'master_publication':
                with plot_tabs[i]:
                    st.image(plot_path, caption=f"{plot_type.replace('_', ' ').title()}", use_container_width=True)

def main():
    """Main application function"""
    initialize_session_state()
    create_header()
    
    # Get sidebar parameters
    params = create_sidebar()
    
    # File upload section
    protein_file, ligand_file, ligand_smiles, ligand_library = create_file_upload_section()
    
    # Main content area
    run_docking = False
    demo_mode = False
    
    if protein_file and (ligand_file or ligand_smiles):
        col1, col2 = st.columns([3, 1])
        
        with col1:
            run_docking = st.button("🚀 Run Docking Analysis", type="primary")
        
        with col2:
            demo_mode = st.button("🎮 Demo Mode", help="Show example results without running actual docking")
    
    elif not protein_file:
        st.info("📁 Please upload a protein receptor file (PDB format) to get started")
        demo_mode = st.button("🎮 Try Demo Mode", help="Show example results without uploading files")
    
    else:
        st.info("💊 Please provide ligand input (file or SMILES string) to proceed")
        demo_mode = st.button("🎮 Try Demo Mode", help="Show example results without uploading files")
    
    if run_docking:
        with st.spinner("Running molecular docking analysis..."):
            try:
                # Create temporary directory for this session
                temp_dir = tempfile.mkdtemp(prefix="pandadock_streamlit_")
                
                # Save uploaded files
                protein_path = os.path.join(temp_dir, protein_file.name)
                with open(protein_path, 'wb') as f:
                    f.write(protein_file.getbuffer())
                
                if ligand_file:
                    ligand_path = os.path.join(temp_dir, ligand_file.name)
                    with open(ligand_path, 'wb') as f:
                        f.write(ligand_file.getbuffer())
                elif ligand_smiles:
                    # Create temporary SDF file from SMILES
                    ligand_path = os.path.join(temp_dir, "ligand.smi")
                    with open(ligand_path, 'w') as f:
                        f.write(ligand_smiles)
                
                # Create configuration
                config = create_config_from_params(params, protein_file, ligand_file, ligand_smiles, temp_dir)
                
                # Run docking
                results, generated_files = run_docking_analysis(config, protein_path, ligand_path)
                
                if results:
                    st.session_state.docking_results = results
                    st.session_state.generated_files = generated_files or {}
                    st.session_state.output_dir = temp_dir
                    
                    st.success("🎉 Docking analysis completed successfully!")
                else:
                    st.error("❌ Docking analysis failed. Please check your input files and parameters.")
                    
            except Exception as e:
                st.error(f"❌ Error during docking analysis: {str(e)}")
                st.info("💡 Try using Demo Mode to see the interface functionality")
    
    elif demo_mode:
        st.info("🎮 Demo Mode: Showing example results")
        # Create demo results
        demo_results = create_demo_results()
        st.session_state.docking_results = demo_results
        st.session_state.generated_files = {}
        st.session_state.output_dir = tempfile.gettempdir()
        st.success("🎉 Demo results loaded successfully!")
    
    # Display results if available
    if hasattr(st.session_state, 'docking_results') and st.session_state.docking_results:
        generated_files = getattr(st.session_state, 'generated_files', {})
        output_dir = getattr(st.session_state, 'output_dir', tempfile.gettempdir())
        
        create_results_visualization(st.session_state.docking_results, generated_files)
        display_generated_plots(generated_files)
        create_file_downloads(generated_files, output_dir)
    
    # Footer
    st.markdown("---")
    st.markdown("**🐼 PandaDock** - Next-Generation Molecular Docking Platform | Built with Streamlit")
    
    # About section
    with st.expander("ℹ️ About PandaDock"):
        st.markdown("""
        **PandaDock** is a comprehensive molecular docking platform featuring three novel algorithms:
        
        - **PandaCore**: Robust baseline algorithm with evolutionary optimization
        - **PandaML**: Advanced machine learning-based algorithm with superior affinity prediction
        - **PandaPhysics**: Physics-based algorithm specialized for metal coordination and complex chemistry
        
        **Key Features:**
        - Publication-quality visualizations
        - Comprehensive binding metrics (IC50, EC50, binding affinity)
        - Interactive 3D molecular visualization
        - Professional interaction analysis
        - Multi-format output (HTML, CSV, JSON, plots)
        
        **Citation:**
        ```
        @software{pandadock2025,
          title={PandaDock: Next-Generation Molecular Docking with Novel PandaDock Algorithms},
          author={Pritam Kumar Panda},
          year={2025},
          url={https://github.com/pritampanda15/pandadock}
        }
        ```
        """)

if __name__ == "__main__":
    main()