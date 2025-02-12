# UKSH Analysis Workflow

This repository contains a reproducible analysis workflow for PDAC patient data. The project is containerized using Docker, ensuring a consistent environment for running the analysis.

## Project Structure

- **Dockerfile:** Contains instructions to build the Docker image.
- **requirements.txt:** Lists all Python dependencies.
- **main.py:** The main script that executes the workflow.
- **functions.py:** Helper functions used in the analysis.
- **Output:** Contains already generated outputs.
- **main_flow.py** Python notebook cointaing code and pre-generated results. 
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

### 3. Output
After executing, output will be generated in terminal and output files will be saved in your current directory.