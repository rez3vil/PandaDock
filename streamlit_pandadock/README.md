# PandaDock Streamlit Web Interface

A comprehensive web interface for PandaDock molecular docking built with Streamlit.

## Features

### 🎯 Complete Docking Workflow
- **Three PandaDock Algorithms**: PandaCore, PandaML, and PandaPhysics
- **Flexible Input Options**: PDB files, SMILES strings, ligand libraries
- **Advanced Configuration**: Grid box settings, flexible docking, scoring parameters
- **Real-time Progress**: Live updates during docking analysis

### 📊 Advanced Visualization
- **Interactive Plots**: Score distributions, binding affinity analysis, IC50 correlations
- **Publication-Quality Figures**: Master publication plots, binding metrics analysis
- **3D Molecular Visualization**: Interactive 3D pose viewers
- **Professional Reports**: HTML, CSV, and JSON export formats

### 🔧 User-Friendly Interface
- **Intuitive Sidebar**: Easy parameter configuration
- **Drag-and-Drop Upload**: Simple file handling
- **Real-time Feedback**: Progress indicators and status updates
- **Comprehensive Downloads**: ZIP archives with all results

## Installation

### Prerequisites
```bash
cd /Users/pritam/PandaDock
pip install -e .
```

### Install Streamlit Dependencies
```bash
cd streamlit_pandadock
pip install -r requirements.txt
```

## Usage

### Start the Application
```bash
cd streamlit_pandadock
streamlit run app.py
```

The application will open in your default web browser at `http://localhost:8501`.

### Basic Workflow

1. **Configure Parameters** (Sidebar)
   - Choose docking algorithm (Balanced/Precise/Fast)
   - Set scoring function (PandaML/PandaPhysics/PandaCore)
   - Adjust docking parameters (poses, exhaustiveness, etc.)
   - Configure grid box settings
   - Set flexible docking options

2. **Upload Files**
   - Upload protein receptor (PDB format)
   - Choose ligand input:
     - Single ligand file (SDF, MOL2, PDB)
     - SMILES string
     - Ligand library file

3. **Run Analysis**
   - Click "🚀 Run Docking Analysis"
   - Monitor progress in real-time
   - View results automatically

4. **Analyze Results**
   - Interactive visualization tabs
   - Statistical analysis plots
   - Binding metrics dashboard
   - Results data table

5. **Download Results**
   - Individual files (plots, reports)
   - Complete ZIP archive
   - Multiple format options

## Interface Sections

### 🔧 Configuration Sidebar
- **Docking Algorithm**: Select from three PandaDock modes
- **Docking Parameters**: Poses, exhaustiveness, energy range
- **Grid Box Settings**: Center coordinates and box size
- **Flexible Docking**: Residue flexibility options
- **Advanced Options**: ML rescoring, random seed

### 📁 Input Files Section
- **Protein Upload**: PDB receptor files
- **Ligand Options**: Files, SMILES, or libraries
- **File Validation**: Automatic format checking
- **Status Indicators**: Upload confirmation

### 📊 Results Analysis
- **Summary Metrics**: Quick overview cards
- **Interactive Tabs**:
  - Score Distribution: Histogram and scatter plots
  - Binding Affinity: Bar charts and correlations
  - IC50 Analysis: Potency distribution
  - Results Table: Sortable data grid

### 🎨 Visualization Gallery
- **Master Publication Plot**: Comprehensive dashboard
- **Individual Analysis Plots**: Specialized visualizations
- **Interactive Display**: Zoom, pan, and explore
- **High-Resolution Images**: Publication-ready quality

### 📥 Download Center
- **ZIP Archive**: All results in one package
- **Individual Files**: Selective downloads
- **Multiple Formats**: PNG, PDF, HTML, CSV, JSON
- **Organized Structure**: Logical file organization

## Advanced Features

### Multi-Algorithm Support
```python
# Algorithms automatically selected based on mode:
# - Balanced → PandaML (ML-based scoring)
# - Precise → PandaPhysics (Physics-based)
# - Fast → PandaCore (Evolutionary optimization)
```

### Comprehensive Output Generation
- HTML reports with interactive 3D viewers
- CSV data files for further analysis
- JSON results for programmatic access
- High-resolution publication plots
- Professional interaction maps

### Real-time Analysis
- Progress tracking during docking
- Live parameter validation
- Immediate result visualization
- Interactive plot updates

## Output Files

### Generated Results
- `pandadock_results.csv` - Tabular results data
- `pandadock_report.html` - Interactive HTML report
- `pandadock_report.json` - Structured JSON data
- `master_publication.png` - Comprehensive analysis plot
- `binding_metrics_analysis.png` - Binding analysis
- `ic50_ec50_analysis.png` - Potency analysis
- `score_distribution_analysis.png` - Statistical analysis

### Pose Files
- `poses/` directory containing:
  - Individual pose PDB files
  - Ligand SDF files
  - Complex structures
  - Summary CSV file

### Visualization Files
- Professional 2D interaction maps
- Interactive 3D HTML visualizations
- Statistical analysis plots
- Publication-quality figures

## Configuration Options

### Docking Modes
- **Balanced**: Best overall performance with ML scoring
- **Precise**: Detailed physics-based analysis
- **Fast**: High-throughput evolutionary optimization

### Scoring Functions
- **PandaML**: Machine learning-based scoring
- **PandaPhysics**: Physics-based energy calculation
- **PandaCore**: Empirical scoring functions

### Grid Box Settings
- **Auto-detection**: Automatic binding site identification
- **Manual**: Custom center coordinates and box size
- **Flexible sizing**: Adjustable dimensions

## Technical Details

### Performance
- **Parallel Processing**: Multi-threading support
- **Memory Efficient**: Optimized for large datasets
- **Scalable**: Handles multiple ligands
- **Responsive**: Real-time interface updates

### File Format Support
- **Protein**: PDB
- **Ligand**: SDF, MOL2, PDB, SMILES
- **Library**: SMI, TXT, CSV
- **Output**: HTML, CSV, JSON, PNG, PDF

### Browser Compatibility
- Chrome/Chromium (recommended)
- Firefox
- Safari
- Edge

## Troubleshooting

### Common Issues
1. **File Upload Errors**: Check file format and size
2. **Docking Failures**: Verify protein/ligand compatibility
3. **Visualization Issues**: Ensure browser JavaScript is enabled
4. **Performance**: Close other applications for better performance

### Error Messages
- Clear error descriptions in interface
- Detailed logging for debugging
- Helpful suggestions for resolution

## Examples

### Basic Docking
1. Upload protein receptor (PDB)
2. Enter ligand SMILES: `CCO`
3. Select "Balanced" mode
4. Click "Run Docking Analysis"
5. Download results

### Advanced Analysis
1. Configure flexible residues: `HIS57,SER195`
2. Set custom grid box coordinates
3. Enable ML rescoring
4. Generate comprehensive plots
5. Export all formats

## Citation

```bibtex
@software{pandadock2025,
  title={PandaDock: Next-Generation Molecular Docking with Novel PandaDock Algorithms},
  author={Pritam Kumar Panda},
  year={2025},
  url={https://github.com/pritampanda15/pandadock}
}
```

## License

MIT License - see the main PandaDock repository for details.

## Support

- **Documentation**: [PandaDock ReadTheDocs](https://pandadock.readthedocs.io/)
- **Issues**: [GitHub Issues](https://github.com/pritampanda15/pandadock/issues)
- **Email**: pritam@stanford.edu

---

**🐼 Dock Smarter. Discover Faster. Visualize Better.**