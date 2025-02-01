# 🧪 SMILES2Mol

**Convert SMILES notation to molecular structures with serverless AWS backend and GitHub Pages frontend.**

[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-blue?logo=github)](https://hyperverseglobalconsulting.github.io/smiles2mol/webpage.html)
[![AWS API Gateway](https://img.shields.io/badge/Endpoint-API%20Gateway-orange)](https://projects.vizeet.me/smiles2mol)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

![Architecture Diagram](smiles2mol-arch-diagram.png) <!-- Add your architecture diagram -->

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

## 🌐 Architecture

### 1. Static Web Hosting
```mermaid
graph LR 
  A[GitHub Repository] --> B[GitHub Pages]
  B --> C[User's Browser]
```
### 2. Backend Processing Pipeline
```mermaid
graph TD
  A[User] -->|SMILES Input| B[GitHub Pages]
  B -->|API Request| C[AWS API Gateway]
  C -->|Trigger| D[Lambda Function]
  D -->|RDKit Processing| E[Return JPG]
  E -->|Display| B
```
### 3. Domain & Certificate Management
1. **ACM Certificate**: Request SSL certificate for `projects.vizeet.me` (us-east-2 region)
2. **API Gateway**: Configure custom domain with certificate
3. **Route53**: Create CNAME record pointing to API Gateway DNS

```mermaid
sequenceDiagram 
  participant A as ACM
  participant B as API Gateway
  participant C as Route53
  A->>A: Request Certificate
  A->>B: Associate with API Gateway
  B->>C: Create CNAME Record
  C->>B: Domain Validation
```

## 🛠️ Implementation

### Frontend (Static Hosting)
```bash
# Clone repository
git clone https://github.com/hyperverseglobalconsulting/smiles2mol.git
```
# Enable GitHub Pages:
1. Go to Repository Settings → Pages
2. Select branch (main/master) and root folder
3. Enforce HTTPS
