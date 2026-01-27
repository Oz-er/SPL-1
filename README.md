# Cross-Project Software Defect Prediction via Transfer Learning

A high-performance C++ implementation of **Transfer Component Analysis (TCA)** for predicting software defects across different projects. This tool leverages domain adaptation to transfer knowledge from a labeled source project (e.g., Apache) to an unlabeled target project (e.g., Zxing/Safe), overcoming the problem of distributional shift in software metrics.

![C++](https://img.shields.io/badge/Language-C++-00599C?style=flat-square)
![Type](https://img.shields.io/badge/Type-From%20Scratch%20Implementation-green?style=flat-square)
![Domain](https://img.shields.io/badge/Domain-Machine%20Learning%20%2F%20SE-blue?style=flat-square)

## 📌 Project Overview

In software engineering, historical defect data is often unavailable for new projects (Target Domain). Traditional machine learning fails here because the data distribution of the new project differs from the old projects (Source Domain).

This project implements **Transfer Learning** algorithms to align these distributions mathematically.

**Key Algorithms Implemented:**
* **TCA (Transfer Component Analysis):** Minimizes the Maximum Mean Discrepancy (MMD) between domains in a latent kernel space.
* **BDA (Balanced Distribution Adaptation):** *(Planned)* Extends TCA to balance marginal and conditional distributions.
* **CORAL (Correlation Alignment):** *(Planned)* Aligns the second-order statistics (covariance) of the source and target.
