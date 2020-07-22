# nmr_autoprocessing
Automated peak alignment for 1H-NMR metabolomics using machine learning.


## computing environment setup
For portability, we will use conda environments. We will also use Jupyter notebooks for exploratory coding prior to converting code to scripts, functions, Classes, etc., so a Jupyter kernel will need to be associated with the conda environment.

First, clone the repository. When you use `git clone`, a new folder that has the same name as the repository is generated in the current directory that you are navigated to within your terminal/commandline:

```
> git clone https://github.com/medlocklab/nmr_autoprocessing.git
```

Create the conda environment; we'll use a conda environment `yaml` file, which specifies the name of the environment we're creating and all dependencies:

```
> conda env create environment.yaml
```

activate the conda environment; you'll need to do this whenever you want to use or modify the environment (e.g., install new packages):

```
> conda activate nmr_autoprocessing
```

The `environment.yaml` configuration file includes the `nb_conda` package, which improves access of conda environments from Jupyter; after activating the environment, you should be able to start notebooks containing a kernel for the `nmr_autoprocessing` environment from within Jupyter. Start JupyterLab to check this and run/modify any notebooks:

```
> jupyter lab
```

