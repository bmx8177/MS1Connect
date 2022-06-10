# MS1Connect

**9 June 2022: Please note that this version of MS1Connect is not fully functional at this time. We hope to upload the remaining code within a few month. If you want me to try MS1Connect on your data please contact me and I will be happy to help.**

MS1Connect is a tool that scores the similarity between a pair of mass spectrometry runs. This task is particularly challenging because data can be acquired under different experimental procotols. MS1Connect solves this problem by framing the problem as a maximum bipartite matching problem and by only using data fom intact peptide (MS1) scans.

We are currently developing the documentation for MS1Connect. Feel free to reach out with any questions in the meantime.

## Citing
If you use MS1Connect in your work please cite:
>https://www.biorxiv.org/content/10.1101/2022.01.12.476125v1

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
