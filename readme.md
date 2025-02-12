# UKSH Analysis Workflow

This repository contains a reproducible analysis workflow for PDAC patient data.

## Project Structure
- **Output:** Contains already generated outputs.
- **Dockerfile:** Contains instructions to build the Docker image.
- **Methods_Resuls:** Contains a report detailing the methods used and results acheived. 
- **functions.py:** Helper functions used in the analysis.
- **main.py:** The main script that executes the workflow.
- **main_flow.py** Python notebook cointaing code and pre-generated results.
- **requirements.txt:** Lists all Python dependencies.
## How to Run the Analysis

### 1. Download the Docker Image

Pull the pre-built image from docker Hub:

```bash
docker pull mmofty/uksh_analysis_workflow:latest
```

### 2. Execute the Docker Container
Run the container with the following command:

```bash
docker run --rm -v "$(pwd):/app" mmofty/uksh_analysis_workflow:latest
```
Data sets are already incuded in the docker image. 


### 3. Output
After executing, output will be generated in terminal and output files will be saved in your current directory.
