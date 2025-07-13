# Setting Up ReadTheDocs for PandaDock

This guide will help you set up ReadTheDocs hosting for the PandaDock documentation.

## Prerequisites

1. **GitHub Repository**: Your PandaDock code should be in a public GitHub repository
2. **ReadTheDocs Account**: Free account at [readthedocs.org](https://readthedocs.org)
3. **Documentation Files**: Sphinx documentation (already configured in `docs/` directory)

## Step-by-Step Setup

### 1. Create ReadTheDocs Account

1. Go to [readthedocs.org](https://readthedocs.org)
2. Click "Sign up" in the top right
3. Choose "Sign up with GitHub" for easy integration
4. Authorize ReadTheDocs to access your GitHub repositories

### 2. Import Your Project

1. After logging in, click "Import a Project"
2. You should see your GitHub repositories listed
3. Find "PandaDock" and click the "+" button to import
4. If you don't see it, click "Import Manually" and enter:
   - **Repository URL**: `https://github.com/your-username/PandaDock`
   - **Project Name**: `PandaDock`
   - **Description**: `Modular, Multi-Strategy, High-Performance Molecular Docking Software`

### 3. Configure Project Settings

1. Go to your project dashboard
2. Click "Admin" → "Settings"
3. Configure these settings:

#### Basic Settings:
- **Name**: PandaDock
- **Description**: Modular, Multi-Strategy, High-Performance Molecular Docking Software
- **Repository URL**: https://github.com/your-username/PandaDock
- **Default branch**: main
- **Default version**: latest

#### Advanced Settings:
- **Configuration file**: `.readthedocs.yaml` (already created)
- **Python version**: 3.11
- **Documentation type**: Sphinx Html
- **Sphinx configuration file**: `docs/source/conf.py`

### 4. Build Settings

The project is configured to build automatically using the `.readthedocs.yaml` file which:
- Uses Python 3.11
- Installs dependencies from `docs/requirements.txt`
- Builds HTML, PDF, and EPUB formats
- Includes benchmark images and analysis

### 5. Custom Domain (Optional)

If you want a custom domain like `docs.pandadock.org`:

1. Go to Admin → Domains
2. Add your custom domain
3. Set up CNAME record: `docs.pandadock.org` → `pandadock.readthedocs.io`

### 6. Webhook Setup (Automatic)

ReadTheDocs automatically creates webhooks when you import via GitHub, so documentation rebuilds automatically when you push changes.

## Testing the Build

### Local Testing

Before pushing, test documentation locally:

```bash
cd docs/
make html
open _build/html/index.html  # macOS
# or
xdg-open _build/html/index.html  # Linux
```

### ReadTheDocs Build

1. Go to your project dashboard
2. Click "Build version" to trigger a manual build
3. Monitor the build logs for any errors
4. Once successful, click "View docs" to see your live documentation

## Current Documentation Structure

The documentation includes:

### User Guide:
- Installation and quickstart
- Configuration and docking modes
- Metal docking and analysis
- Comprehensive benchmarks

### Tutorials:
- Basic docking tutorial
- ML-enhanced docking (PandaML)
- Metal docking tutorial
- Comprehensive benchmark tutorial
- Virtual screening

### API Reference:
- Complete API documentation
- Code examples
- Function references

### Benchmarks:
- Comprehensive results on 5,316 complexes
- Metal vs non-metal analysis (1,982 vs 3,334 complexes)
- Performance comparisons
- Statistical analysis

## Customization Options

### Theme Customization

The current Sphinx theme can be customized in `docs/source/conf.py`:

```python
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'canonical_url': '',
    'analytics_id': 'UA-XXXXXXX-1',  # Add Google Analytics
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    'style_nav_header_background': 'white',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}
```

### Adding Google Analytics

1. Get a Google Analytics tracking ID
2. Add to `docs/source/conf.py`:
   ```python
   html_theme_options = {
       'analytics_id': 'UA-XXXXXXX-1',
   }
   ```

### Custom CSS/JS

Add custom styling in `docs/source/_static/custom.css` and reference in `conf.py`:

```python
html_static_path = ['_static']
html_css_files = [
    'custom.css',
]
```

## Troubleshooting

### Common Issues:

#### 1. Build Fails - Missing Dependencies
**Error**: `ModuleNotFoundError`
**Solution**: Update `docs/requirements.txt` with missing packages

#### 2. Images Not Loading
**Error**: Benchmark plots not displaying
**Solution**: Ensure image paths are correct and files are committed to Git

#### 3. Sphinx Warnings
**Error**: Various Sphinx warnings
**Solution**: Check `docs/source/conf.py` configuration and RST syntax

#### 4. PDF Build Fails
**Error**: LaTeX compilation issues
**Solution**: Disable PDF builds in `.readthedocs.yaml` if not needed:
```yaml
formats: []  # Remove PDF/EPUB if causing issues
```

### Getting Help:

1. **ReadTheDocs Support**: [docs.readthedocs.io/en/stable/support.html](https://docs.readthedocs.io/en/stable/support.html)
2. **Sphinx Documentation**: [www.sphinx-doc.org](https://www.sphinx-doc.org)
3. **GitHub Issues**: Create issues in your repository for project-specific problems

## Expected URLs

Once set up, your documentation will be available at:

- **Main site**: `https://pandadock.readthedocs.io/`
- **Latest version**: `https://pandadock.readthedocs.io/en/latest/`
- **Stable version**: `https://pandadock.readthedocs.io/en/stable/`
- **PDF download**: `https://pandadock.readthedocs.io/_/downloads/en/latest/pdf/`

## Maintenance

### Regular Updates:
1. **Keep dependencies updated** in `docs/requirements.txt`
2. **Monitor build status** for any failures
3. **Update documentation** with new features and benchmarks
4. **Review analytics** to understand documentation usage

### Version Management:
- ReadTheDocs automatically builds documentation for different Git branches/tags
- Use semantic versioning for releases (v1.0.0, v1.1.0, etc.)
- Configure which versions to display in ReadTheDocs admin panel

The documentation will automatically rebuild whenever you push changes to GitHub, keeping your docs always up-to-date with the latest PandaDock developments and benchmark results!