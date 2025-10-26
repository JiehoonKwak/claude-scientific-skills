# Scientific Skills API Keys & Configuration Requirements

**Analysis Date**: 2025-10-23
**Total Skills Analyzed**: 88
**Skills Requiring Configuration**: 17 (19%)

---

## Quick Reference: Skills Requiring API Keys

### ✅ MANDATORY Configuration (9 skills)

1. **COSMIC Database** - Academic registration or commercial license
2. **USPTO Database** - API key required
3. **Benchling Integration** - API key or OAuth credentials
4. **DNAnexus Integration** - Account credentials via `dx login`
5. **LabArchives Integration** - API Access Key ID + Password
6. **OMERO Integration** - Server credentials
7. **LatchBio Integration** - Account credentials via `latch login`
8. **biomni** - LLM API keys (Anthropic, OpenAI, etc.)
9. **paper-2-web** - OpenAI API key
10. **pymatgen** - Materials Project API key

### ⚠️ RECOMMENDED/OPTIONAL Configuration (7 skills)

11. **FDA Database** - Optional (increases rate limits)
12. **ClinPGx Database** - Optional (for heavy use)
13. **biopython** - Email required, API key optional
14. **bioservices** - Email for some NCBI services
15. **gget** - API keys for specific modules only

---

## Detailed Breakdown by Category

## 1. Scientific Databases (4/24 require configuration)

### ✅ COSMIC Database
- **Status**: MANDATORY
- **Academic Users**: Free registration
  - URL: https://cancer.sanger.ac.uk/cosmic/register
  - Credentials: Email + password
- **Commercial Users**: License required through QIAGEN
- **Usage**: Download access for mutation data

### ✅ USPTO Database
- **Status**: MANDATORY
- **Registration**: https://account.uspto.gov/api-manager/
- **Environment Variable**: `USPTO_API_KEY`
- **Usage**: All USPTO API endpoints

### ⚠️ FDA Database
- **Status**: RECOMMENDED
- **Registration**: https://open.fda.gov/apis/authentication/
- **Rate Limits**:
  - Without key: 240 req/min, 1,000/day
  - With key: 240 req/min, 120,000/day
- **Usage**: Enhanced rate limits for heavy use

### ⚠️ ClinPGx Database
- **Status**: OPTIONAL
- **Contact**: api@clinpgx.org (for substantial use)
- **Rate Limit**: 2 requests/second maximum
- **Usage**: Basic access works without key

### ❌ No Configuration Required (20 databases)

- AlphaFold Database
- ChEMBL Database
- ClinicalTrials.gov
- ClinVar Database
- ENA Database
- Ensembl Database
- Gene Database (NCBI)
- GEO Database
- GWAS Catalog
- HMDB Database
- KEGG Database (academic use only)
- Open Targets
- Metabolomics Workbench
- PDB Database
- PubChem Database
- PubMed Database
- Reactome Database
- STRING Database
- UniProt Database
- ZINC Database

---

## 2. Scientific Integrations (6/7 require configuration)

### ✅ Benchling Integration
- **Status**: MANDATORY
- **Authentication Methods**:
  1. API Key (from Profile Settings)
  2. OAuth Client Credentials
  3. OIDC (advanced)
- **Required**: Tenant URL (e.g., https://your-tenant.benchling.com)
- **Security**: Store in environment variables

### ✅ DNAnexus Integration
- **Status**: MANDATORY
- **Setup**: `dx login` command
- **Required**: DNAnexus account credentials
- **Verification**: `dx whoami`

### ✅ LabArchives Integration
- **Status**: MANDATORY
- **Required**:
  - API Access Key ID
  - API Access Password
  - User email + external applications password
  - Regional API endpoint (US/Australia/UK)
- **Prerequisites**: Enterprise LabArchives license with API access
- **Configuration File**: `config.yaml`
  ```yaml
  api_url: https://api.labarchives.com/api
  access_key_id: YOUR_ACCESS_KEY_ID
  access_password: YOUR_ACCESS_PASSWORD
  ```

### ✅ OMERO Integration
- **Status**: MANDATORY
- **Required**:
  - OMERO server hostname
  - Port number
  - Username
  - Password
- **Connection**:
  ```python
  conn = BlitzGateway(username, password, host=host, port=port)
  ```

### ✅ LatchBio Integration
- **Status**: MANDATORY
- **Setup**: `latch login` command
- **Prerequisites**:
  - Latch account credentials
  - Docker installed and running

### ❌ Opentrons Integration
- **Status**: NOT REQUIRED (for protocol development)
- **Note**: Physical robot needed for execution, but not for development/simulation

---

## 3. Scientific Packages (7/47 require configuration)

### ✅ biomni
- **Status**: MANDATORY
- **Required**: `ANTHROPIC_API_KEY` or other LLM provider keys:
  - OpenAI
  - Azure
  - Google
  - Groq
  - AWS Bedrock
- **Configuration**: `.env` file or environment variables
- **Optional**: Google API keys for enhanced features

### ⚠️ biopython
- **Status**: EMAIL REQUIRED, API KEY OPTIONAL
- **Required**: Email address for NCBI Entrez API
  ```python
  Entrez.email = "your_email@example.com"
  ```
- **Optional**: NCBI API key for higher rate limits
  ```python
  Entrez.api_key = "your_api_key"
  ```

### ⚠️ bioservices
- **Status**: MOSTLY NOT REQUIRED
- **Note**: Some NCBI services require email address (not a key)

### ⚠️ gget
- **Status**: MODULE-SPECIFIC
- **Requires API Keys**:
  - `gget gpt`: OpenAI API key
  - `gget cosmic`: COSMIC credentials
- **No API Key**: Most other modules work without keys

### ✅ paper-2-web
- **Status**: MANDATORY
- **Required**: `OPENAI_API_KEY` (all operations)
- **Optional**:
  - `GOOGLE_API_KEY`
  - `GOOGLE_CSE_ID` (for logo search)

### ✅ pymatgen
- **Status**: MANDATORY
- **Required**: Materials Project API key
- **Environment Variable**: `MP_API_KEY`
- **Setup**: `export MP_API_KEY="your_key_here"`
- **Registration**: https://next-gen.materialsproject.org/

### ⚠️ cellxgene-census
- **Status**: NOT REQUIRED
- **Note**: Large data downloads (~11GB on first use)

### ⚠️ diffdock
- **Status**: NOT REQUIRED
- **Note**: Downloads ~500MB model parameters on first run, GPU recommended

### ❌ No Configuration Required (40 packages)

- arboreto
- anndata
- astropy
- cobrapy
- dask
- datamol
- deepchem
- deeptools
- etetoolkit
- flowio
- matplotlib
- matchms
- medchem
- molfeat
- polars
- pydeseq2
- pymc
- pymoo
- pyopenms
- pysam
- pytorch-lightning
- pytdc
- rdkit
- reportlab
- scanpy
- scikit-learn
- scikit-bio
- statsmodels
- seaborn
- torch_geometric
- torchdrug
- transformers
- zarr-python
- umap-learn
- tooluniverse

---

## 4. Scientific Thinking Skills (0/12 require configuration)

### ❌ All Skills - No Configuration Required

**Document Skills**:
- docx (requires LibreOffice, pandoc, Node.js packages)
- pdf (requires pypdf, pdfplumber, reportlab)
- pptx (requires LibreOffice, markitdown)
- xlsx (requires LibreOffice, openpyxl)

**Thinking Skills**:
- hypothesis-generation
- exploratory-data-analysis
- peer-review
- scientific-brainstorming
- scientific-critical-thinking
- scientific-visualization
- scientific-writing
- statistical-analysis

**Note**: These skills require software packages (Python libraries, LibreOffice, etc.) but NO API keys or authentication credentials.

---

## Summary Statistics

| Category | Requires Config | No Config | Total | % Requiring Config |
|----------|----------------|-----------|-------|-------------------|
| Databases | 4 | 20 | 24 | 17% |
| Integrations | 6 | 1 | 7 | 86% |
| Packages | 7 | 40 | 47 | 15% |
| Thinking | 0 | 12 | 12 | 0% |
| **TOTAL** | **17** | **71** | **88** | **19%** |

---

## Configuration Priority Levels

### 🔴 High Priority (Commonly Used)
1. **biopython** - NCBI email (widely used in bioinformatics)
2. **pymatgen** - Materials Project key (materials science)
3. **biomni** - LLM keys (AI-enhanced workflows)

### 🟡 Medium Priority (Domain-Specific)
4. **COSMIC** - Cancer research
5. **Benchling** - Lab data management
6. **DNAnexus** - Genomics cloud platform

### 🟢 Low Priority (Specialized Use Cases)
7. **USPTO** - Patent searching
8. **LabArchives** - Electronic lab notebooks
9. **OMERO** - Microscopy image management
10. **LatchBio** - Workflow orchestration

---

## Environment Variables Template

Create a `.env` file with required API keys:

```bash
# LLM APIs (for biomni, paper-2-web, gget gpt)
ANTHROPIC_API_KEY=your_anthropic_key_here
OPENAI_API_KEY=your_openai_key_here

# Scientific Databases
USPTO_API_KEY=your_uspto_key_here
FDA_API_KEY=your_fda_key_here  # Optional
MP_API_KEY=your_materials_project_key_here

# NCBI Services (for biopython, bioservices)
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your_ncbi_key_here  # Optional

# COSMIC Database
COSMIC_EMAIL=your_cosmic_email
COSMIC_PASSWORD=your_cosmic_password

# Integrations
BENCHLING_API_KEY=your_benchling_key_here
BENCHLING_TENANT_URL=https://your-tenant.benchling.com

LABARCHIVES_ACCESS_KEY_ID=your_labarchives_key_id
LABARCHIVES_ACCESS_PASSWORD=your_labarchives_password
LABARCHIVES_API_URL=https://api.labarchives.com/api

OMERO_HOST=your_omero_server
OMERO_PORT=4064
OMERO_USERNAME=your_username
OMERO_PASSWORD=your_password

# Optional - Google Services
GOOGLE_API_KEY=your_google_key_here
GOOGLE_CSE_ID=your_cse_id_here
```

---

## Setup Checklist

### Phase 1: Essential APIs (No Cost)
- [ ] Set NCBI email for biopython
- [ ] Register for Materials Project API key
- [ ] Register for FDA API key (optional but recommended)

### Phase 2: LLM APIs (Paid Services)
- [ ] Set up Anthropic API key for biomni
- [ ] Set up OpenAI API key for paper-2-web, gget gpt

### Phase 3: Specialized Databases (Registration Required)
- [ ] Register for COSMIC (academic) or obtain commercial license
- [ ] Register for USPTO API key
- [ ] Obtain NCBI API key for higher rate limits (optional)

### Phase 4: Platform Integrations (Enterprise/Institutional)
- [ ] Set up Benchling credentials
- [ ] Set up DNAnexus account (`dx login`)
- [ ] Configure LabArchives API access
- [ ] Configure OMERO server connection
- [ ] Set up LatchBio account (`latch login`)

---

## Security Best Practices

1. **Never commit API keys to version control**
   - Add `.env` to `.gitignore`
   - Use environment variables or secure vaults

2. **Use separate keys for development/production**
   - Isolate testing from production workloads

3. **Rotate keys periodically**
   - Especially for high-value APIs (LLM providers, commercial databases)

4. **Monitor API usage**
   - Track rate limits and costs
   - Set up alerts for unusual activity

5. **Use minimal permissions**
   - Request only necessary scopes for OAuth
   - Use read-only keys where possible

---

## Next Steps

1. **Identify required skills** - Review which scientific skills you need for your work
2. **Prioritize API setup** - Start with free/essential APIs before paid services
3. **Create `.env` file** - Use the template above for your configuration
4. **Test connections** - Verify each API key works before full deployment
5. **Document internal processes** - Track which keys are used in which projects

---

## Resources

### Free API Keys
- Materials Project: https://next-gen.materialsproject.org/
- FDA OpenAPI: https://open.fda.gov/apis/authentication/
- NCBI E-utilities: https://www.ncbi.nlm.nih.gov/account/
- COSMIC (academic): https://cancer.sanger.ac.uk/cosmic/register

### Paid/Enterprise APIs
- Anthropic: https://console.anthropic.com/
- OpenAI: https://platform.openai.com/
- Benchling: https://www.benchling.com/
- DNAnexus: https://www.dnanexus.com/
- LabArchives: https://www.labarchives.com/
- LatchBio: https://latch.bio/

### Documentation
- USPTO API: https://developer.uspto.gov/
- NCBI Entrez: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- BioPython Tutorial: https://biopython.org/DIST/docs/tutorial/Tutorial.html

---

*Generated by Claude Code analysis of 88 scientific skills*
