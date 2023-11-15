# <Project Name>

This code accompanies the paper entitled:<Project Name>


To reproduce these results, follow these steps:

1.  Clone this GitHub repo (i.e. the page you are on):
    * ```git clone https://github.com/NIDAP-Community/<Project Name>.git```

2.  The input files for this pipeline will be available upon request. Please reach out to the authors before continue to following steps

3.  Install docker and build the docker container:
    * Navigate to the cloned repository directory. 
    * Move to the ./Docker_file/ directory of this repo

4.  Build the container:
    * ```docker build --tag <Project Name> .```

5.  Navidate to the cloned repository directory, Run the conainer by mounting the ./src/ directory of the repo to /tmp/ in the container:
    * ```docker run -ti -v $(pwd)/src:/mnt <Project Name>```
    
6.  Run the following code.
    * ```cd /mnt```
    * ```bash run_pipeline.sh```

