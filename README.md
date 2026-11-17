# ChemSASA

**Chemical Functional Group Solvent-Accessible Surface Area (ChemSASA)** is a Python-based tool for calculating solvent-accessible surface area (SASA) from protein structures and molecular dynamics (MD) trajectories. In addition to total and atom-wise SASA, ChemSASA resolves SASA by chemical functional group based on the functional-group classification used in the **solute partitioning model (SPM)** [1,2].

ChemSASA uses [MDAnalysis](https://www.mdanalysis.org/) [3,4] for reading molecular structures and trajectories and [FreeSASA](https://freesasa.github.io/) [5] for solvent-accessible surface area calculations. The current ChemSASA workflow uses the Lee-Richards algorithm with a probe radius of 1.4 Å, 100 slices, and the Richards (1977) van der Waals radii [6]. Users can change the probe radius and number of slices or use the Shrake-Rupley algorithm by modifying the corresponding FreeSASA parameters in the `SASA` module. A different set of atomic radii can also be used by providing a different FreeSASA classifier file.

## Installation

### Clone the repository

```bash
git clone <repository-url>
cd ChemSASA
```

### Create the Python environment

The required Python packages and tested package versions are specified in `environment.yaml`.

```bash
conda env create -f environment.yaml
conda activate chemsasa
```

## Parameters

The parameter files used by ChemSASA are provided in `config-files/`.

- **`freesasa_richards1977.config`** contains the atomic radii used by FreeSASA for SASA calculations, based on the van der Waals radii of Richards [6]. Users can provide a different FreeSASA configuration file to use another set of atomic radii.
- **`res-func-grp.json`** contains the chemical functional-group definitions used by chemsasa.py to assign amino acid atoms to their corresponding functional groups. The definitions are based on the functional-group classification used in the SPM [1,2].

## Example notebook and data

`chemsasa.py` contains the complete ChemSASA implementation in a single Python file. The individual classes are also provided as separate modules in `chemsasa_modules/` for use in custom workflows.

- **`res_func_grp`** contains the chemical functional-group definitions used to assign amino acid atoms to their corresponding functional groups.
- **`UniqueResidues`** reads the topology and counts the atoms assigned to each chemical functional group.
- **`SASA`** performs the SASA calculations for structures and MD trajectories and assigns the calculated atomic SASA values to their corresponding chemical functional groups.

`examples/ChemSASA_example.ipynb` provides detailed example workflows for a PDB structure and an MD trajectory. The examples use ubiquitin (`1UBQ`) as the PDB structure and a short energy-minimized 1UBQ trajectory generated for the trajectory calculation.

The main repository organization is:

```text
ChemSASA/
├── chemsasa.py
├── chemsasa_modules/
├── config-files/
├── examples/
│   ├── ChemSASA_example.ipynb
│   └── example-data/
├── environment.yaml
└── README.md
```

The example notebook includes the complete workflow and sanity checks for both calculations. Reference output files are provided with the example data in the corresponding `output-files/` directories. 

## References

1. Guinn, E. J., Pegram, L. M., Capp, M. W., Pollock, M. N., & Record, M. T. Jr. (2011). Quantifying why urea is a protein denaturant, whereas glycine betaine is a protein stabilizer. *Proceedings of the National Academy of Sciences*, **108**(41), 16932–16937. https://doi.org/10.1073/pnas.1109372108

2. Cheng, X., Guinn, E. J., Buechel, E., Wong, R., Sengupta, R., Shkel, I. A., & Record, M. T. Jr. (2016). Basis of protein stabilization by K glutamate: Unfavorable interactions with carbon, oxygen groups. *Biophysical Journal*, **111**(9), 1854–1865. https://doi.org/10.1016/j.bpj.2016.08.050

3. Michaud-Agrawal, N., Denning, E. J., Woolf, T. B., & Beckstein, O. (2011). MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *Journal of Computational Chemistry*, **32**, 2319–2327. https://doi.org/10.1002/jcc.21787

4. Gowers, R. J., Linke, M., Barnoud, J., Reddy, T. J. E., Melo, M. N., Seyler, S. L., Dotson, D. L., Domański, J., Buchoux, S., Kenney, I. M., & Beckstein, O. (2016). MDAnalysis: A Python package for the rapid analysis of molecular dynamics simulations. *Proceedings of the 15th Python in Science Conference*, 98–105. https://doi.org/10.25080/Majora-629e541a-00e

5. Mitternacht, S. (2016). FreeSASA: An open source C library for solvent accessible surface area calculations. *F1000Research*, **5**, 189. https://doi.org/10.12688/f1000research.7931.1

6. Richards, F. M. (1977). Areas, volumes, packing, and protein structure. *Annual Review of Biophysics and Bioengineering*, **6**, 151–176. https://doi.org/10.1146/annurev.bb.06.060177.001055
