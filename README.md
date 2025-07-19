# T2-MRI-Myelin-Preterm-Analysis

This repository contains a project developed for the COMP0118 module, focusing on the analysis of T2 relaxometry MRI data to study brain microstructure and estimate myelin content in preterm-born adolescents.

## 📁 Repository Structure

- **`COMP0118-Analysis of MRI T2 Relaxometry- CODE/`**  
  Contains all the MATLAB scripts for:
  - Single- and multi-compartment T2 model fitting
  - Myelin water fraction (MWF) estimation
  - Signal evaluation and residual analysis
  - Model comparison (AIC/BIC)
  - Region-based analysis using segmentation masks

- **`COMP0118 Presentation T2 Relaxometry.pptx`**  
  A PowerPoint presentation summarizing the goals, methodology, and preliminary results from the first part of the project.

- **`COMP0118-Analysis of MRI T2 Relaxometry.pdf`**  
  The final project report, written in scientific article format. Includes background, methods, results with figures and tables, and discussion.

- **`T2_MRI_coursework.pdf`**  
  The coursework specification document detailing the project objectives, tasks, and dataset description.

## 🔍 Project Summary

The project explores quantitative T2 relaxometry for brain MRI, aiming to characterize white matter microstructure and estimate tissue-specific T2 values. It implements:

- Mono- and multi-exponential signal decay modeling
- Model fitting via least squares (LS, NLLS, NNLS)
- Multi-compartment modeling (2–4 compartments, fixed and free T2)
- Myelin-specific signal extraction (MWF)
- Model evaluation using residuals and information criteria

## 🧠 Clinical Context

This analysis is part of a larger research effort to assess long-term effects of extreme prematurity on brain development. Estimating myelin content through T2 relaxometry provides insight into white matter integrity and potential links to cognitive outcomes.

## 💻 Requirements

- MATLAB R2021a or later
- NIfTI I/O functions (e.g., `load_nii`)
- Optimization Toolbox (for `fmincon`)

## 📜 License

This project is for academic use only, distributed under an educational fair-use policy.

## Author

Francesco Seracini  
MSc Computer Science - Artificial Intelligence  
Politecnico di Milano - University College London  

