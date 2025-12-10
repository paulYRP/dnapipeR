--------------
# DNAm EPIC Array Analysis Workflow
--------------

This repository contains a reproducible pipeline for DNA methylation (DNAm) analysis using Illumina MethylationEPIC v2.0 arrays. It includes preprocessing, quality control, phenotype merging, and statistical modeling (GLM and LMM) for any phenotype using both `minfi`, `watermelon` and `ewastools` frameworks. The pipeline is modular, HPC-compatible, and fully parameterized via command-line arguments. 

--------------
## Articles:
- [**A Pilot Epigenome-Wide Study of Posttraumatic Growth: Identifying Novel Candidates for Future Research**](https://www.mdpi.com/2075-4655/9/4/39)

--------------
## Tutorials:
- [**DNA Methylation Tutorial**](https://paulYRP.github.io/2025-cpgpneurogenomics-workshop/tutorial.html)
- [**Getting Started**](https://github.com/paulYRP/DNAm_ArrayWorkflow/wiki/Getting-Started)
- [**Requirements**](https://github.com/paulYRP/DNAm_ArrayWorkflow/wiki/Requirements)

--------------
## News

### **10/12/2025**
A first functional version of **dnapipeR** has been created in the `package` branch.
* **GitHub installation support** enabled via:
  ```r
  devtools::install_github("paulYRP/dnapipeR@package")
  ```
  
### 28/10/2025
- Added support for **multiple models** using the `MODEL` and `MODELS` variables inside the Makefile.
  - `MODEL ?= model1`, sets the default model if none is specified.
  - `MODELS = model1 model2 model3`, defines all models to run in parallel.
- Each model runs in **isolated directories**. 
- Parallel execution added using:
  - `make models`, runs the **entire pipeline** for all models in parallel.
  - `make f3_models`, runs **Steps 1–3** in parallel for all models.
  - `make f4_models`, runs **Steps 1–4** in parallel for all models.
  - `make f3lme_models`, runs **Steps 1–3 + LME** in parallel for all models.    
- Individual model runs remain supported:
  - `make all MODEL=model1`, runs the full pipeline for a single model.
  - `make f3 MODEL=model1`, runs the first 3 steps for that model only.

### 21/08/2025
- It is compatible with all types of tissues using [ewastools](https://hhhh5.github.io/ewastools/articles/exemplary_ewas.html) 
- `ctrlsva` added from [ENmix](https://www.bioconductor.org/packages/devel/bioc/vignettes/ENmix/inst/doc/ENmix.html) 
- `adjusted_funnorm` added from [wateRmelon](https://www.bioconductor.org/packages/release/bioc/vignettes/wateRmelon/inst/doc/wateRmelon.html)
- `make -j1 f3` added to execute only the first three steps. Output:
  - Metrics (beta, M)
  - Objects (RGSet, MSet ...)
  - CSV file with cell estimation
  - Surrogate variable report
  - CSV file and ZIP file compatible with [ClockFundation](https://dnamage.clockfoundation.org/)
- `make -j1 f4` added to execute only the first 4 steps. Output:
  - annotatedGLM.csv
- `make -j1 f3lme` added to execute only the first 3 + LMER. Output:
  - annotatedLME.csv
- `DNA.pdf`, automic report generated with a summary of all steps.
- Interaction added to `methylationGLM_T1`.
- `methylationGLM_T1` and `methylationGLMM_T1T2` extract categorical and numerical coefficients to the annotatedGLM.csv and annotatedLME.csv. 
