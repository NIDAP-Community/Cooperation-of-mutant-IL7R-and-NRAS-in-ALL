# Cooperation-of-mutant-IL7R-and-NRAS-in-ALL

This code accompanies the paper entitled:

This is a repository for the code and results from this manuscript: Winer H, Wenqing L, Rodrigues G, Gower T, Meyer TJ, Hixon J, Durum SK. Mechanism of co-operation of mutant IL-7Ralpha and mutant NRAS in acute lymphoblastic leukemia: role of MYC. Haematologica. (Accepted, November 2023).

To reproduce these results, follow these steps:

1.  Clone this GitHub repo (i.e. the page you are on):
    * ```git clone https://github.com/NIDAP-Community/Cooperation-of-mutant-IL7R-and-NRAS-in-ALL```

2.  The input files for this pipeline will be available upon request. Please reach out to the authors before continue to following steps

3.  Install [docker](https://docs.docker.com/get-docker/) and build the docker container:
    * Navigate to the cloned repository directory. 
    * Move to the ./Docker_file/ directory of this repo

4.  Build the container:
    * ```docker build --tag cooperation-of-mutant-il7r-and-nras-in-all .```

5.  Navidate to the cloned repository directory, Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container:
    * ```docker run -ti -v $(pwd)/src:/mnt cooperation-of-mutant-il7r-and-nras-in-all```
    
6.  Run the following code.
    * ```cd /mnt```
    * ```bash run_pipeline.sh```

