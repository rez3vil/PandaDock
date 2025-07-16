# Contributing to PandaDock

Thank you for your interest in contributing to PandaDock! This document provides guidelines and information for contributors.

## Table of Contents

1. [Code of Conduct](#code-of-conduct)
2. [Getting Started](#getting-started)
3. [Development Setup](#development-setup)
4. [Contributing Guidelines](#contributing-guidelines)
5. [Pull Request Process](#pull-request-process)
6. [Issue Guidelines](#issue-guidelines)
7. [Development Workflow](#development-workflow)
8. [Testing](#testing)
9. [Documentation](#documentation)
10. [Code Style](#code-style)
11. [Algorithm Contributions](#algorithm-contributions)

## Code of Conduct

This project and everyone participating in it is governed by our commitment to creating a welcoming and inclusive environment. By participating, you are expected to uphold professional standards of conduct.

### Our Standards

- **Be respectful**: Treat all contributors with respect and courtesy
- **Be inclusive**: Welcome newcomers and help them get started
- **Be collaborative**: Work together constructively and share knowledge
- **Be patient**: Understand that people have different skill levels and backgrounds
- **Be constructive**: Provide helpful feedback and suggestions

## Getting Started

### Prerequisites

- Python 3.8 or higher
- Git for version control
- Basic understanding of molecular docking concepts
- Familiarity with scientific Python libraries (NumPy, SciPy, etc.)

### Quick Start

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/YOUR_USERNAME/PandaDock.git
   cd PandaDock
   ```
3. **Set up development environment** (see [Development Setup](#development-setup))
4. **Create a feature branch** for your contribution
5. **Make your changes** following our guidelines
6. **Submit a pull request** when ready

## Development Setup

### Environment Setup

1. **Create a virtual environment**:
   ```bash
   python -m venv pandadock-dev
   source pandadock-dev/bin/activate  # On Windows: pandadock-dev\Scripts\activate
   ```

2. **Install in development mode**:
   ```bash
   pip install -e .[dev,test,docs]
   ```

3. **Install pre-commit hooks**:
   ```bash
   pre-commit install
   ```

### Verify Installation

```bash
# Run basic tests
pytest tests/test_basic.py

# Check code style
black --check pandadock/
flake8 pandadock/

# Verify imports
python -c "import pandadock; print(pandadock.__version__)"
```

## Contributing Guidelines

### Types of Contributions

We welcome various types of contributions:

- **🐛 Bug fixes**: Fix issues and improve stability
- **✨ New features**: Add new functionality or algorithms
- **📚 Documentation**: Improve docs, examples, and tutorials
- **🔧 Performance**: Optimize existing code
- **🧪 Tests**: Add or improve test coverage
- **🎨 Visualization**: Enhance plotting and analysis capabilities
- **🔬 Algorithms**: Contribute new docking or scoring algorithms

### Contribution Process

1. **Check existing issues** to see if your idea is already being worked on
2. **Open an issue** to discuss significant changes before starting work
3. **Follow the development workflow** outlined below
4. **Write tests** for new functionality
5. **Update documentation** as needed
6. **Submit a pull request** following our template

## Pull Request Process

### Before Submitting

- [ ] Code follows project style guidelines
- [ ] All tests pass locally
- [ ] New tests added for new functionality
- [ ] Documentation updated if needed
- [ ] CHANGELOG.md updated for significant changes
- [ ] Pre-commit hooks pass

### PR Template

```markdown
## Description
Brief description of changes and motivation.

## Type of Change
- [ ] Bug fix (non-breaking change fixing an issue)
- [ ] New feature (non-breaking change adding functionality)
- [ ] Breaking change (fix or feature causing existing functionality to change)
- [ ] Documentation update

## Testing
- [ ] Existing tests pass
- [ ] New tests added
- [ ] Manual testing performed

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] Tests added/updated
```

### Review Process

1. **Automated checks** must pass (CI/CD pipeline)
2. **Code review** by maintainers
3. **Testing** on different platforms if needed
4. **Documentation review** for user-facing changes
5. **Approval** and merge by maintainers

## Issue Guidelines

### Bug Reports

Use the bug report template and include:

- **Clear description** of the issue
- **Steps to reproduce** the problem
- **Expected vs actual behavior**
- **Environment details** (OS, Python version, dependencies)
- **Error messages** and stack traces
- **Minimal example** demonstrating the issue

### Feature Requests

Use the feature request template and include:

- **Clear description** of the proposed feature
- **Use case** and motivation
- **Proposed implementation** approach
- **Alternatives considered**
- **Impact assessment** on existing functionality

### Algorithm Proposals

For new docking or scoring algorithms:

- **Scientific background** and references
- **Performance characteristics** and benchmarks
- **Integration approach** with existing framework
- **Implementation plan** and timeline

## Development Workflow

### Branch Naming

- `feature/algorithm-name` - New algorithms or major features
- `fix/issue-description` - Bug fixes
- `docs/topic` - Documentation improvements
- `test/component` - Test additions/improvements
- `refactor/component` - Code refactoring

### Commit Messages

Follow conventional commit format:

```
type(scope): description

[optional body]

[optional footer]
```

Types: `feat`, `fix`, `docs`, `test`, `refactor`, `perf`, `chore`

Examples:
```
feat(algorithms): add PandaML diffusion sampling
fix(scoring): correct VdW energy calculation
docs(api): update algorithm selection guide
test(physics): add metal coordination tests
```

### Git Workflow

1. **Create feature branch**:
   ```bash
   git checkout -b feature/my-feature
   ```

2. **Make commits** with clear messages:
   ```bash
   git add .
   git commit -m "feat(component): description"
   ```

3. **Keep branch updated**:
   ```bash
   git fetch upstream
   git rebase upstream/main
   ```

4. **Push and create PR**:
   ```bash
   git push origin feature/my-feature
   ```

## Testing

### Test Structure

```
tests/
├── test_algorithms.py      # Algorithm testing
├── test_scoring.py         # Scoring function tests
├── test_io.py             # Input/output tests
├── test_integration.py    # Integration tests
└── data/                  # Test data files
```

### Running Tests

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_algorithms.py

# Run with coverage
pytest --cov=pandadock

# Run integration tests
pytest -m integration

# Skip slow tests
pytest -m "not slow"
```

### Writing Tests

- **Use descriptive names**: `test_pandaml_generates_valid_poses`
- **Follow AAA pattern**: Arrange, Act, Assert
- **Test edge cases**: Empty inputs, invalid parameters
- **Use fixtures**: For common test data and setup
- **Mock external dependencies**: Use `pytest-mock` for external services

### Test Example

```python
def test_pandaml_confidence_prediction():
    """Test PandaML confidence prediction functionality."""
    # Arrange
    config = PandaDockConfig()
    engine = MLEngine(config)
    
    # Act
    confidence = engine.predict_confidence(protein_features, ligand_features, pose)
    
    # Assert
    assert 0.0 <= confidence <= 1.0
    assert isinstance(confidence, float)
```

## Documentation

### Documentation Types

- **API Documentation**: Docstrings following NumPy style
- **User Guides**: Tutorials and how-to guides
- **Algorithm Documentation**: Technical implementation details
- **Examples**: Jupyter notebooks and Python scripts

### Docstring Style

Follow NumPy docstring conventions:

```python
def calculate_energy(ligand_coords, protein_coords=None):
    """
    Calculate total binding energy for a ligand pose.
    
    Parameters
    ----------
    ligand_coords : np.ndarray
        Ligand coordinates with shape (n_atoms, 3)
    protein_coords : np.ndarray, optional
        Protein coordinates with shape (m_atoms, 3)
        
    Returns
    -------
    float
        Total binding energy in kcal/mol (negative for favorable binding)
        
    Examples
    --------
    >>> coords = np.array([[0, 0, 0], [1, 1, 1]])
    >>> energy = calculate_energy(coords)
    >>> energy < 0  # Favorable binding
    True
    """
```

### Building Documentation

```bash
# Install documentation dependencies
pip install -e .[docs]

# Build HTML documentation
cd docs/
make html

# Preview documentation
python -m http.server 8000 -d _build/html/
```

## Code Style

### Python Style

- **PEP 8** compliance with 88-character line length
- **Type hints** for all public functions
- **Black** for code formatting
- **Ruff** for linting
- **isort** for import sorting

### Pre-commit Hooks

Automatic code formatting and linting:

```bash
# Install pre-commit
pre-commit install

# Run manually
pre-commit run --all-files
```

### Configuration

See `pyproject.toml` for tool configurations:
- Black formatting
- Ruff linting rules
- isort import sorting
- MyPy type checking

## Algorithm Contributions

### Adding New Algorithms

1. **Inherit from DockingEngine**:
   ```python
   class MyAlgorithm(DockingEngine):
       def dock(self, protein_file, ligand_file):
           # Implementation
   ```

2. **Follow naming conventions**:
   - File: `pandadock/docking/my_algorithm.py`
   - Class: `MyAlgorithm`
   - Config: `my_algorithm` scoring function

3. **Implement required methods**:
   - `dock()`: Main docking method
   - `score()`: Pose scoring
   - `get_engine_info()`: Algorithm metadata

4. **Add comprehensive tests**:
   - Unit tests for core functionality
   - Integration tests with scoring functions
   - Performance benchmarks

5. **Update documentation**:
   - Algorithm wiki entry
   - Usage examples
   - Performance characteristics

### Scoring Functions

1. **Add to ScoringFunctions class**:
   ```python
   def calculate_my_energy(self, ligand_coords, protein_coords=None):
       """Calculate energy using my method."""
   ```

2. **Update energy dispatcher**:
   ```python
   elif scoring_function == 'my_algorithm':
       return self.calculate_my_energy(ligand_coords, protein_coords)
   ```

3. **Add configuration support**:
   - Update config schema
   - Add algorithm-specific parameters
   - Document configuration options

## Getting Help

### Resources

- **Documentation**: https://pandadock.readthedocs.io/
- **Algorithm Wiki**: `PANDADOCK_ALGORITHMS_WIKI.md`
- **Examples**: `examples/` directory
- **Issues**: GitHub issue tracker

### Communication

- **GitHub Issues**: For bugs and feature requests
- **GitHub Discussions**: For questions and ideas
- **Email**: For sensitive issues (pritam@stanford.edu)

### Maintainer Response

- **Bug reports**: Within 48 hours
- **Feature requests**: Within 1 week
- **Pull requests**: Within 1 week
- **Documentation**: Within 72 hours

## Recognition

Contributors are recognized in:

- **CHANGELOG.md**: For significant contributions
- **README.md**: In the acknowledgments section
- **Documentation**: In contributor guides
- **Releases**: In release notes

## License

By contributing to PandaDock, you agree that your contributions will be licensed under the MIT License.

---

## Quick Reference

### Essential Commands

```bash
# Setup development environment
python -m venv pandadock-dev
source pandadock-dev/bin/activate
pip install -e .[dev,test,docs]
pre-commit install

# Code quality checks
black pandadock/
ruff pandadock/
mypy pandadock/

# Testing
pytest
pytest --cov=pandadock

# Documentation
cd docs/ && make html
```

### File Structure

```
pandadock/
├── __init__.py              # Package initialization
├── config.py                # Configuration system
├── docking/                 # Docking algorithms
│   ├── ga_engine.py        # PandaCore algorithm
│   ├── ml_engine.py        # PandaML algorithm
│   └── physics_engine.py   # PandaPhysics algorithm
├── scoring/                 # Scoring functions
├── io/                      # Input/output handling
├── reports/                 # Analysis and reporting
└── utils/                   # Utility functions
```

Thank you for contributing to PandaDock! 🐼