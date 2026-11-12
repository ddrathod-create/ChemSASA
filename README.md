# ChemSASA

ChemSASA is a Python-based tool for calculating solvent-accessible surface area (SASA) from protein structures and molecular dynamics (MD) trajectories. In addition to total and atom-wise SASA, ChemSASA resolves solvent-accessible surface area by chemical functional group based on the functional-group classification used in the **solute partitioning model (SPM)** [1,2].

ChemSASA uses [MDAnalysis](https://www.mdanalysis.org/) [3,4] for handling molecular structures and trajectories and [FreeSASA](https://freesasa.github.io/) [5] for solvent-accessible surface area calculations.

## Installation

### Clone the repository

```bash
git clone <repository-url>
cd ChemSASA
```

### Create the Python environment

The required Python packages and tested package versions are specified in `environment.yaml`. Create the Conda environment using:

```bash
conda env create -f environment.yaml
```

Activate the environment:

```bash
conda activate chemsasa
```

## Parameters

The `params/` directory contains the parameter files used by ChemSASA for assigning atoms to chemical functional groups and defining atomic radii for SASA calculations.

- **`res-func-grp.json`** defines the assignment of protein atoms to the chemical functional groups used by ChemSASA. The functional-group classification is based on the solute partitioning model (SPM) [1,2].

- **`freesasa_richards1977.config`** contains the atomic radii used by FreeSASA for SASA calculations, based on the van der Waals radii of Richards [6].

## Usage

An example workflow using the individual ChemSASA modules is provided in `ChemSASA_example.ipynb`. The notebook demonstrates how to set up and run ChemSASA for both a single protein structure and an MD trajectory.

The repository is organized as follows:

```text
ChemSASA/
├── chemsasa.py
├── ChemSASA_example.ipynb
├── chemsasa-modules/
│   └── ...
├── params/
│   ├── res-func-grp.json
│   └── freesasa_richards1977.config
├── environment.yaml
└── README.md
```

The individual ChemSASA classes are provided in `chemsasa-modules/` and are used in the example notebook. The complete implementation is also available in `chemsasa.py`.

## References

1. Guinn, E. J., Pegram, L. M., Capp, M. W., Pollock, M. N., & Record, M. T. Jr. (2011). Quantifying why urea is a protein denaturant, whereas glycine betaine is a protein stabilizer. *Proceedings of the National Academy of Sciences*, **108**(41), 16932–16937. https://doi.org/10.1073/pnas.1109372108

2. Cheng, X., Guinn, E. J., Buechel, E., Wong, R., Sengupta, R., Shkel, I. A., & Record, M. T. Jr. (2016). Basis of protein stabilization by K glutamate: Unfavorable interactions with carbon, oxygen groups. *Biophysical Journal*, **111**(9), 1854–1865. https://doi.org/10.1016/j.bpj.2016.08.050

3. Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., & Beckstein, O. (2011). MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *Journal of Computational Chemistry*, **32**, 2319–2327. https://doi.org/10.1002/jcc.21787

4. Gowers, R. J., Linke, M., Barnoud, J., Reddy, T. J. E., Melo, M. N., Seyler, S. L., Dotson, D. L., Domański, J., Buchoux, S., Kenney, I. M., & Beckstein, O. (2016). MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. *Proceedings of the 15th Python in Science Conference*, 98–105. https://doi.org/10.25080/Majora-629e541a-00e

5. Mitternacht, S. (2016). FreeSASA: An open source C library for solvent accessible surface area calculations. *F1000Research*, **5**, 189. https://doi.org/10.12688/f1000research.7931.1

6. Richards, F. M. (1977). Areas, volumes, packing, and protein structure. *Annual Review of Biophysics and Bioengineering*, **6**, 151–176. https://doi.org/10.1146/annurev.bb.06.060177.001055
