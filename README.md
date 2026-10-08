# ChemSASA

ChemSASA is a Python-based tool for calculating solvent-accessible surface area (SASA) from protein structures and molecular dynamics (MD) trajectories. In addition to total and atom-wise SASA, SASAmap resolves solvent-accessible surface area by chemical functional group.

ChemSASA uses [MDAnalysis](https://www.mdanalysis.org/) for handling molecular structures and trajectories and [FreeSASA](https://freesasa.github.io/) for solvent-accessible surface area calculations.

## Installation

### Clone the repository

```bash
git clone <repository-url>
cd ChemSASA
```

### Create the Python environment

The required Python packages and tested package versions are specified in `environment.yaml`. Create a Conda environment from this file using:

```bash
conda env create -f environment.yaml
```

Activate the environment:

```bash
conda activate chemsasa-env
```
