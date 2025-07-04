# Replication Folder for "Effects of Monetary Policy on the Wealth Inequality"

This repository contains the replication materials for the paper currently under Revise & Resubmit at *Economics Letters*.  
Please note that this is a preliminary version and may be subject to updates.

(Last updated: July 4, 2025)  
Yeongwoong Do

---

## 1. Data Description

The `raw data` folder includes the following:

- **scf**: Survey of Consumer Finances (SCF) data  
- **dfa-networth-levels-detail.xlsx**: Distributional Financial Accounts (DFA) data  
- **FRED_data**: Macroeconomic variables (including interest rates) manually collected from FRED  
- **monetary-policy-surprises-data.xlsx**: Monetary policy shocks from Bauer and Swanson (2023)  
- **Wealth_Realtimeinequality**: Real-Time Inequality (RTI) data used in Appendix analysis related to weak instrument robustness  

---

## 2. Replication Instructions

Run the `STATA do` files and `MATLAB m` files sequentially following the `Step_xx` numbering.  
The process will automatically generate the following folders:

- `Figure`: Figures corresponding to the main text and appendix  
- `STATA_dta`: Intermediate data files  
- `Table`: Tables included in the main text and appendix  

Below is a summary of key files directly linked to the main text results:

| Output                                                          | Step File                | Software |
|-----------------------------------------------------------------|--------------------------|----------|
| **Figure 1**: Portfolio composition by asset class              | `Step_02.do`             | STATA    |
| **Table 1**: Wealth Gini decomposition by asset components      | `Step_03.do`             | STATA    |
| **Figure 2**: Impulse responses of wealth Gini to monetary tightening | `Step_07.m`           | MATLAB   |
| **Figure 3**: Impulse responses of component Gini and shares to monetary tightening | `Step_10.do` | STATA    |
| **Table 2**: IRF decomposition                                  | `Step_11.do`             | STATA    |

---

## Notes

- Some raw data files (e.g., SCF microdata) may require separate access or user agreements depending on their source.
- All codes are provided for academic use and transparency.


