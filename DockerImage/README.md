# Griffin Docker
Griffin docker container aims to streamline the setup process for Griffin and enhance the user experience by encapsulating all dependencies and required configurations. 

## Getting Started
To get started with the Griffin Docker container, please follow the instructions below:
1. Change your working direcotry
```
cd Griffin/DockerImage
```
2. Build the Docker image `griffin-docker:v0.2.0` using Dockerfile
```
docker build --platform linux/amd64 -t griffin-docker:v0.2.0 .
```
3. Run the Griffin Docker container
```
docker run -it --rm griffin-docker:v0.2.0
```

## Run Griffin Demo
If you would like to run the [Griffin Demo](https://github.com/adoebley/Griffin/wiki), it is recommended to use the 
Dockerfile located in `Griffin/DockerImage/demo_DockerImage`. This Dockerfile contains all necessary reference files 
and a bash script that helps you get started with running the Griffin Demo.

1. Change your working direcotry
```
cd Griffin/DockerImage/demo_DockerImage
```
2. Build the Docker image `griffin-docker:demo` using demo Dockerfile
```
docker build --platform linux/amd64 -t griffin-docker:demo .
```
3. Run the Demo Griffin Docker container
```
docker run -it --rm griffin-docker:demo
```
4. Run Griffin Demo
Once inside the container, you can start using Griffin by following the Griffin wiki guide.
The demo user guide has been summarized into `run-demo.sh`. The bash script will dry run 
GC correction and griffin_nucleosome_profiling snakefile. 
```
bash run-demo.sh
```

## Happy nucleosome profiling!
