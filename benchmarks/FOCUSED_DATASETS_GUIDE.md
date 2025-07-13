# 🎯 Focused Molecular Docking Datasets Guide

A comprehensive guide to downloading and using high-quality, focused datasets for superior binding affinity prediction performance.

---

## 🏆 Tier 1: Premium Datasets (Highest Quality)

### **1. CASF (Comparative Assessment of Scoring Functions)**
**Best for**: Overall scoring function validation  
**Expected R²**: 0.4-0.6  
**Size**: ~285 high-quality complexes

#### Download & Setup:
```bash
# Download CASF-2016 (most recent)
wget http://www.pdbbind.org.cn/download/CASF-2016.tar.gz
tar -xzf CASF-2016.tar.gz

# Structure:
# CASF-2016/
#   ├── coreset/           # 285 complexes
#   ├── power_docking/     # For docking power assessment
#   ├── power_ranking/     # For ranking power assessment
#   └── power_screening/   # For screening power assessment
```

**Official Website**: http://www.pdbbind.org.cn/casf.php  
**Paper**: J. Chem. Inf. Model. 2019, 59, 895-913

---

### **2. PDBbind Refined Set**
**Best for**: High-quality diverse validation  
**Expected R²**: 0.3-0.5  
**Size**: ~5,000 complexes (high resolution + good binding data)

#### Download & Setup:
```bash
# Download PDBbind v2020 (requires registration)
# Visit: http://www.pdbbind.org.cn/download.php

# Register for academic license (free)
# Download: PDBbind_v2020_refined.tar.gz

tar -xzf PDBbind_v2020_refined.tar.gz

# Structure:
# refined-set/
#   ├── index/
#   │   └── INDEX_refined_data.2020
#   ├── 1a30/
#   │   ├── 1a30_protein.pdb
#   │   ├── 1a30_ligand.mol2
#   │   └── 1a30_ligand.sdf
#   └── [4000+ other complexes]
```

**Registration**: http://www.pdbbind.org.cn/download.php  
**License**: Free for academic use

---

## 🧬 Tier 2: Protein Family-Specific Datasets

### **3. Kinase Inhibitor Dataset (ChEMBL)**
**Best for**: Kinase drug discovery  
**Expected R²**: 0.5-0.7  
**Size**: ~2,000-5,000 complexes

#### Download Script:
```python
# save as download_kinase_data.py
from chembl_webresource_client.new_client import new_client

# Get kinase targets
target = new_client.target
kinase_targets = target.filter(
    target_type="PROTEIN FAMILY",
    pref_name__icontains="kinase"
).only(['target_chembl_id', 'pref_name'])

# Get kinase activities
activity = new_client.activity
kinase_activities = activity.filter(
    target_chembl_id__in=[t['target_chembl_id'] for t in kinase_targets[:100]],
    standard_type__in=['Ki', 'Kd', 'IC50'],
    standard_value__isnull=False,
    pchembl_value__isnull=False
)

# Save to CSV
import pandas as pd
df = pd.DataFrame(kinase_activities)
df.to_csv('kinase_dataset.csv', index=False)
```

#### Run:
```bash
pip install chembl_webresource_client
python download_kinase_data.py
```

**Alternative**: Use preprocessed datasets from:  
- **BindingDB**: https://www.bindingdb.org/bind/chemsearch/marvin/SDFdownload.jsp
- **ChEMBL Interface**: https://www.ebi.ac.uk/chembl/

---

### **4. GPCR Ligand Dataset**
**Best for**: GPCR drug discovery  
**Expected R²**: 0.4-0.6  
**Size**: ~1,000-3,000 complexes

#### Download from GPCRdb:
```bash
# Visit GPCRdb structure browser
# https://gpcrdb.org/structure/

# Download all GPCR structures with ligands
wget "https://gpcrdb.org/structure/download_structures" \
     -O gpcr_structures.zip

# Or use their API
curl -X GET "https://gpcrdb.org/services/structure/" \
     -H "accept: application/json" > gpcr_list.json
```

**Official Sources**:
- **GPCRdb**: https://gpcrdb.org/
- **GPCR-EXP**: http://zhanglab.ccmb.med.umich.edu/GPCR-EXP/
- **GLIDA**: https://glida.fudan.edu.cn/

---

### **5. Protease Inhibitor Dataset**
**Best for**: Protease drug discovery  
**Expected R²**: 0.5-0.7  
**Size**: ~500-1,500 complexes

#### Manual Collection:
```bash
# Download from MEROPS (protease database)
# https://www.ebi.ac.uk/merops/

# Search PDB for protease structures
curl "https://search.rcsb.org/rcsbsearch/v2/query" \
     -H "Content-Type: application/json" \
     -d '{
       "query": {
         "type": "terminal",
         "service": "text",
         "parameters": {
           "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
           "operator": "contains_words",
           "value": "protease"
         }
       },
       "return_type": "entry",
       "request_options": {
         "return_all_hits": true
       }
     }'
```

**Curated Sources**:
- **HIV Protease Database**: https://hivdb.stanford.edu/
- **MEROPS**: https://www.ebi.ac.uk/merops/
- **Protein Data Bank**: https://www.rcsb.org/ (search "protease inhibitor")

---

## 🧪 Tier 3: Specialized & Fragment Datasets

### **6. Astex Diverse Set (Fragment-like)**
**Best for**: Fragment-based drug discovery  
**Expected R²**: 0.6-0.8  
**Size**: 85 high-quality complexes

#### Download:
```bash
# Download from Cambridge Crystallographic Data Centre
# https://www.ccdc.cam.ac.uk/Community/blog/astex-diverse-set/

# Or from supplementary data of original paper
wget https://pubs.acs.org/doi/suppl/10.1021/jm070596t/suppl_file/jm070596t_si_001.pdf

# PDB codes included in the set:
# 1a9u, 1b9o, 1bzc, 1c1b, 1c1u, 1c4u, 1c83, 1c84, 1c87, 1c88
# [... full list in supplementary material]
```

**Paper**: J. Med. Chem. 2007, 50, 5076-5084  
**Focus**: MW < 300 Da, drug-like fragments

---

### **7. DUD-E (Database of Useful Decoys Enhanced)**
**Best for**: Virtual screening validation  
**Expected R²**: 0.4-0.6  
**Size**: ~22,000 actives + decoys across 102 targets

#### Download:
```bash
# Download all targets
wget -r -np -nH --cut-dirs=1 \
     http://dude.docking.org/targets/

# Or download specific targets:
# Kinases
wget http://dude.docking.org/targets/abl1/abl1.tar.gz
wget http://dude.docking.org/targets/cdk2/cdk2.tar.gz

# GPCRs  
wget http://dude.docking.org/targets/adra1a/adra1a.tar.gz
wget http://dude.docking.org/targets/adrb1/adrb1.tar.gz

# Extract
for file in *.tar.gz; do tar -xzf "$file"; done
```

**Website**: http://dude.docking.org/  
**Focus**: Benchmarking virtual screening performance

---

### **8. MOAD (Mother of All Databases)**
**Best for**: Diverse binding site analysis  
**Expected R²**: 0.3-0.5  
**Size**: ~35,000 binding sites

#### Download:
```bash
# Download from University of Michigan
wget http://www.bindingmoad.org/files/biln/every_part_a.zip
wget http://www.bindingmoad.org/files/biln/every_part_b.zip

unzip every_part_a.zip
unzip every_part_b.zip

# Or use their search interface
# http://www.bindingmoad.org/Home/search
```

**Website**: http://www.bindingmoad.org/  
**Focus**: All binding sites in PDB with biological relevance

---

## 🛠️ Dataset Preparation Scripts

### **Universal Dataset Converter**
```python
# save as prepare_dataset.py
import os
import pandas as pd
from pathlib import Path

def prepare_pdbbind_format(input_dir, output_dir):
    """Convert any dataset to PDBbind-like format"""
    
    # Create output structure
    os.makedirs(f"{output_dir}/index", exist_ok=True)
    
    # Process each complex
    index_data = []
    
    for pdb_dir in Path(input_dir).iterdir():
        if pdb_dir.is_dir():
            pdb_code = pdb_dir.name
            
            # Create output directory
            out_dir = Path(output_dir) / pdb_code
            os.makedirs(out_dir, exist_ok=True)
            
            # Copy/convert files
            for file in pdb_dir.iterdir():
                if file.suffix == '.pdb':
                    # Split into protein and ligand if needed
                    split_protein_ligand(file, out_dir, pdb_code)
                elif file.suffix in ['.mol2', '.sdf']:
                    # Copy ligand file
                    shutil.copy(file, out_dir / f"{pdb_code}_ligand{file.suffix}")
            
            # Add to index
            index_data.append({
                'pdb_code': pdb_code,
                'resolution': 'N/A',
                'year': 'N/A', 
                'affinity': 'N/A',  # Will need manual annotation
                'reference': 'N/A',
                'ligand_name': 'UNK'
            })
    
    # Save index file
    df = pd.DataFrame(index_data)
    df.to_csv(f"{output_dir}/index/INDEX_prepared_data.csv", index=False)

def split_protein_ligand(pdb_file, output_dir, pdb_code):
    """Split PDB file into protein and ligand components"""
    protein_lines = []
    ligand_lines = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if line[17:20].strip() in ['HOH', 'WAT']:  # Skip water
                    continue
                elif line.startswith('HETATM'):
                    ligand_lines.append(line)
                else:
                    protein_lines.append(line)
            elif line.startswith(('HEADER', 'CRYST1', 'REMARK')):
                protein_lines.append(line)
    
    # Save protein
    with open(output_dir / f"{pdb_code}_protein.pdb", 'w') as f:
        f.writelines(protein_lines)
        f.write("END\n")
    
    # Save ligand  
    if ligand_lines:
        with open(output_dir / f"{pdb_code}_ligand.pdb", 'w') as f:
            f.writelines(ligand_lines)
            f.write("END\n")

# Usage
if __name__ == "__main__":
    prepare_pdbbind_format("raw_dataset/", "prepared_dataset/")
```

### **Binding Affinity Annotator**
```python
# save as annotate_affinities.py
import pandas as pd
import requests
import time

def get_binding_data_from_chembl(pdb_code):
    """Get binding affinity data from ChEMBL API"""
    
    url = f"https://www.ebi.ac.uk/chembl/api/data/activity"
    params = {
        'target_chembl_id': pdb_code,
        'standard_type__in': 'Ki,Kd,IC50',
        'format': 'json'
    }
    
    try:
        response = requests.get(url, params=params)
        if response.status_code == 200:
            data = response.json()
            if data['activities']:
                # Return best available affinity
                for activity in data['activities']:
                    if activity.get('pchembl_value'):
                        return float(activity['pchembl_value'])
        return None
    except:
        return None

def annotate_dataset(dataset_dir):
    """Add binding affinity annotations to dataset"""
    
    index_file = f"{dataset_dir}/index/INDEX_prepared_data.csv"
    df = pd.read_csv(index_file)
    
    for idx, row in df.iterrows():
        pdb_code = row['pdb_code']
        
        # Try to get binding data
        affinity = get_binding_data_from_chembl(pdb_code)
        
        if affinity:
            df.at[idx, 'affinity'] = affinity
            print(f"✅ {pdb_code}: {affinity} pKd")
        else:
            print(f"❌ {pdb_code}: No binding data found")
        
        time.sleep(0.1)  # Be nice to APIs
    
    # Save updated index
    df.to_csv(index_file, index=False)
    print(f"\n📊 Annotated {len(df[df['affinity'] != 'N/A'])} complexes with binding data")

# Usage
if __name__ == "__main__":
    annotate_dataset("prepared_dataset/")
```

---

## 🚀 Quick Start Guide

### **Option 1: Start with CASF (Recommended)**
```bash
# 1. Download CASF-2016
wget http://www.pdbbind.org.cn/download/CASF-2016.tar.gz
tar -xzf CASF-2016.tar.gz

# 2. Run PandaDock benchmark
cd /Users/pritam/PandaDock
python benchmarks/scripts/comprehensive_benchmark.py \
    --pdbbind_dir CASF-2016/coreset \
    --output_dir casf_benchmark \
    --max_complexes 285

# Expected result: R² = 0.4-0.6
```

### **Option 2: Kinase-Focused Dataset**
```bash
# 1. Download kinase subset from PDBbind
# (Manual curation of kinase complexes from refined set)

# 2. Create kinase-specific directory structure
mkdir kinase_dataset

# 3. Run focused benchmark
python benchmarks/scripts/comprehensive_benchmark.py \
    --pdbbind_dir kinase_dataset \
    --output_dir kinase_benchmark \
    --max_complexes 500

# Expected result: R² = 0.5-0.7
```

### **Option 3: Fragment-like Dataset (Astex)**
```bash
# 1. Get Astex Diverse Set PDB codes
# 2. Download structures using fetch_pdb.py script
# 3. Run benchmark on small, high-quality set

# Expected result: R² = 0.6-0.8
```

---

## 📊 Expected Performance by Dataset

| Dataset | Size | Expected R² | Best For |
|---------|------|-------------|----------|
| **CASF-2016** | 285 | 0.4-0.6 | Overall validation |
| **PDBbind Refined** | 5,000 | 0.3-0.5 | Large-scale robust |
| **Kinase ChEMBL** | 2,000+ | 0.5-0.7 | Drug discovery |
| **GPCR Focused** | 1,000+ | 0.4-0.6 | Membrane proteins |
| **Astex Diverse** | 85 | 0.6-0.8 | Fragment SBDD |
| **DUD-E** | 22,000 | 0.4-0.6 | Virtual screening |

---

## 🔗 Essential Resources

### **Primary Sources**:
- **PDBbind**: http://www.pdbbind.org.cn/
- **ChEMBL**: https://www.ebi.ac.uk/chembl/
- **BindingDB**: https://www.bindingdb.org/
- **Protein Data Bank**: https://www.rcsb.org/

### **Specialized Databases**:
- **GPCRdb**: https://gpcrdb.org/
- **KLIFS (Kinases)**: https://klifs.net/
- **CASTp (Binding Sites)**: http://sts.bioe.uic.edu/castp/
- **ProCognis**: https://www.procognis.com/

### **Academic Papers**:
- **CASF**: J. Chem. Inf. Model. 2019, 59, 895-913
- **PDBbind**: Nucleic Acids Res. 2005, 33, D233-D237  
- **DUD-E**: J. Med. Chem. 2012, 55, 6582-6594
- **Astex**: J. Med. Chem. 2007, 50, 5076-5084

---

## ⚡ Pro Tips

1. **Start Small**: Begin with CASF-2016 (285 complexes) for quick validation
2. **Quality over Quantity**: 500 high-quality complexes > 5,000 mixed quality
3. **Document Everything**: Keep detailed logs of data sources and preprocessing
4. **Version Control**: Track dataset versions and preprocessing steps
5. **Validate Results**: Compare with published benchmarks from literature

---

**Ready to achieve R² > 0.5?** Start with CASF-2016 - it's the gold standard for scoring function validation! 🎯