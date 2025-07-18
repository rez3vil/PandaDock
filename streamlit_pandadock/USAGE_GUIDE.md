# PandaDock Streamlit App - Usage Guide

## 🚀 Quick Start

### 1. Launch the Application
```bash
cd /Users/pritam/PandaDock/streamlit_pandadock
streamlit run app.py
```

Or use the provided script:
```bash
./run.sh
```

### 2. Access the Web Interface
Open your browser and go to: `http://localhost:8501`

## 🎯 Using the Interface

### Step 1: Configure Parameters (Sidebar)
- **Docking Mode**: Choose from Balanced, Precise, or Fast
- **Scoring Function**: Select PandaML, PandaPhysics, or PandaCore
- **Docking Parameters**: Set number of poses, exhaustiveness, energy range
- **Grid Box**: Use auto-detection or set custom coordinates
- **Flexible Docking**: Specify flexible residues
- **Advanced Options**: Enable ML rescoring, set random seed

### Step 2: Upload Files
- **Protein**: Upload PDB file of your receptor
- **Ligand**: Choose from:
  - Single ligand file (SDF, MOL2, PDB)
  - SMILES string (e.g., "CCO" for ethanol)
  - Ligand library file

### Step 3: Run Analysis
- **Run Docking**: Perform actual molecular docking
- **Demo Mode**: See example results without uploading files

### Step 4: View Results
- **Summary Metrics**: Quick overview of results
- **Interactive Plots**: Score distributions, binding affinity, IC50 analysis
- **Results Table**: Detailed data in tabular format
- **Generated Plots**: Publication-quality visualizations

### Step 5: Download Results
- **ZIP Archive**: All results in one package
- **Individual Files**: Select specific plots or reports
- **Multiple Formats**: PNG, PDF, HTML, CSV, JSON

## 🎮 Demo Mode

If you want to test the interface without uploading files:

1. Click "🎮 Try Demo Mode" button
2. View example results with realistic molecular docking data
3. Explore all visualization and download features
4. Perfect for learning the interface or demonstrations

## 📊 Available Visualizations

### Interactive Plots
- **Score Distribution**: Histogram and scatter plots of docking scores
- **Binding Affinity**: Bar charts and correlation analysis
- **IC50 Analysis**: Potency distribution and relationships
- **Results Table**: Sortable data grid with all metrics

### Generated Plots (when available)
- **Master Publication Plot**: Comprehensive analysis dashboard
- **Binding Metrics Analysis**: Detailed binding analysis
- **IC50/EC50 Analysis**: Drug potency analysis
- **Score Distribution Analysis**: Statistical validation

## 🔧 Configuration Options

### Docking Modes
- **Balanced** (PandaML): Best overall performance with ML scoring
- **Precise** (PandaPhysics): Detailed physics-based analysis
- **Fast** (PandaCore): High-throughput evolutionary optimization

### Scoring Functions
- **PandaML**: Machine learning-based scoring
- **PandaPhysics**: Physics-based energy calculation
- **PandaCore**: Empirical scoring functions

### Parameters
- **Number of Poses**: 1-50 poses (default: 10)
- **Exhaustiveness**: 1-32 (default: 8, higher = more thorough)
- **Energy Range**: 1-10 kcal/mol (default: 3.0)
- **Grid Box Size**: 10-50 Å (default: 22.5 Å)

## 📁 File Formats

### Input Files
- **Protein**: PDB format
- **Ligand**: SDF, MOL2, PDB formats
- **SMILES**: Text string (e.g., "CCO", "C1=CC=CC=C1")
- **Library**: SMI, TXT, CSV formats

### Output Files
- **Plots**: PNG, PDF formats
- **Reports**: HTML, CSV, JSON formats
- **Structures**: PDB, SDF formats
- **Archive**: ZIP file with all results

## 🧪 Example Workflows

### Basic Docking
1. Upload protein PDB file
2. Enter ligand SMILES: `CCO`
3. Select "Balanced" mode
4. Click "Run Docking Analysis"
5. View results and download

### Advanced Analysis
1. Upload protein and ligand files
2. Set flexible residues: `HIS57,SER195`
3. Configure custom grid box
4. Enable ML rescoring
5. Generate comprehensive plots
6. Export all formats

### Virtual Screening
1. Upload protein receptor
2. Upload ligand library file
3. Select "Fast" mode for screening
4. Set higher number of poses
5. Analyze top hits

## 🛠️ Troubleshooting

### Common Issues
- **File Upload Errors**: Check file format and size
- **Docking Failures**: Verify protein/ligand compatibility
- **Visualization Issues**: Ensure browser JavaScript is enabled
- **Performance**: Close other applications for better performance

### Error Messages
- Clear descriptions appear in the interface
- Use Demo Mode to test interface functionality
- Check console for detailed error logs

### Tips
- Use Demo Mode to familiarize yourself with the interface
- Start with default parameters for initial testing
- Gradually adjust parameters based on your needs
- Save successful configurations for future use

## 🎨 Interface Features

### Modern Design
- Clean, professional layout
- Intuitive navigation
- Responsive design for all screen sizes
- Custom color scheme matching PandaDock branding

### Interactive Elements
- Real-time parameter validation
- Progress indicators during analysis
- Hover tooltips for guidance
- Expandable help sections

### Data Visualization
- Interactive Plotly charts
- Zoomable and pannable plots
- Downloadable visualizations
- Professional color schemes

## 🔄 Updates and Maintenance

### Regular Updates
- New visualization features
- Enhanced error handling
- Performance improvements
- Additional file format support

### Feedback
- Report issues via GitHub
- Suggest new features
- Share usage examples
- Contribute improvements

## 📞 Support

- **GitHub Issues**: Report bugs and request features
- **Email**: pritam@stanford.edu
- **Documentation**: Comprehensive guides and tutorials
- **Community**: User discussions and examples

---

**🐼 Happy Docking with PandaDock Streamlit!**