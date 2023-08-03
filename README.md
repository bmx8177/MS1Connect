# MS1Connect

MS1Connect is a tool that scores the similarity between a pair of mass
spectrometry runs. This task is particularly challenging because data can be
acquired under different experimental procotols. MS1Connect solves this problem
by framing the problem as a maximum bipartite matching problem and by only using
data fom intact peptide (MS1) scans.

We are currently developing the documentation for MS1Connect. Feel free to reach
out with any questions in the meantime.

## Install
Before you run ms1connect.py you must run the makefile. You can do this by
running the following command.
```
make
```

## Running MS1Connect
Before running MS1Connect we highly suggest you create a new enviroment (such as
conda).

The inputs to MS1Connect is a set of MS1 features files. MS1Connect can generate
these files for you using pyOpenMS by providing a folder of mzML files. If you
prefer to use your own MS1 feature detection method you can instead provide a
folder of MS1 feature files that contains the following columns: m/z, intensity,
retention time, noramlized retention time (proportion of TIC), and charge.

Running MS1Connect requires running the ms1connect.py script. You can see how to
run this script by running the following command.
```
python ms1connect.py -h
```

## Citing
If you use MS1Connect in your work please cite:
>Lin A, Deatherage Kaiser BL, Hutchison JR, Bilmes JA, Noble WS. MS1Connect: a
>mass spectrometry run similarity measure. Bioinformatics. 2023 Feb 3;39(2)

The manuscript describing MS1Connect can be found <a
href="https://academic.oup.com/bioinformatics/article/39/2/btad058/7005198">here</a>.


## Dependencies
MS1Connect requires the following:
- Python
- C++ (gcc)
- Singularity/Docker

In addition the following Python packages are required:
- pyOpenMS
- Pandas
- Numba
- Scipy
- Matplotlib
- Scikit-learn
- Seaborn 
