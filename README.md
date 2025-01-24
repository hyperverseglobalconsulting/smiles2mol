# 🧪 SMILES2Mol

**Convert SMILES notation to molecular structures with serverless AWS backend and GitHub Pages frontend.**

[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-blue?logo=github)](https://hyperverseglobalconsulting.github.io/smiles2mol/)
[![AWS API Gateway](https://img.shields.io/badge/Endpoint-API%20Gateway-orange)](https://projects.vizeet.me/smiles2mol)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

![Architecture Diagram]([architecture.png](smiles2mol-arch-diagram.png) <!-- Add your architecture diagram -->

---

## 📖 Overview

A cloud-native tool to convert SMILES strings into 2D molecular structures:
- **Frontend**: Static web interface hosted on GitHub Pages
- **Backend**: Serverless AWS infrastructure for SMILES processing

---

## 🌐 Architecture

### Frontend
- **Hosting**: GitHub Pages (zero-cost static hosting)
- **Tech Stack**:
  - HTML/JavaScript for UI
  - Simple form interface with image display
  - Fetches data from AWS API Gateway
- **URL**: [https://hyperverseglobalconsulting.github.io/smiles2mol/webpage.html](https://hyperverseglobalconsulting.github.io/smiles2mol/webpage.html)

### Backend
- **API Gateway**: 
  - Custom domain: `projects.vizeet.me/smiles2mol`
  - HTTPS enabled via AWS ACM-managed SSL certificate
- **Lambda Function**:
  - Python-based SMILES processing
  - Uses RDKit for molecular structure generation
- **Infrastructure**:
  - Deployed with AWS CDK/CloudFormation
  - Auto-scaling and pay-per-use pricing

```mermaid
graph TD
  A[User] -->|SMILES Input| B[GitHub Pages]
  B -->|API Request| C[AWS API Gateway]
  C -->|Trigger| D[Lambda Function]
  D -->|RDKit Processing| E[Return SVG/PNG]
  E -->|Display| B
