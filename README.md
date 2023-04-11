# Oncology BIomaRker Discovery (OncoBird)

Prerequisites
--
The OncoBird shiny application runs in a local virtual [docker](https://docs.docker.com/) environment. A Bioconductor package is currently in progress, but is already publicly available, for instructions regarding its installation see `code/OncoBird/`.

Running the OncoBird shiny application with docker
--
After starting docker, you can pull the docker image for OncoBird by typing:
```
docker pull aljoshoh/oncobird_shiny
```
Finally, run the image by:
```
docker run --rm -p 3838:3838 aljoshoh/oncobird_shiny
```
Now you can navigate to `127.0.0.1:3838` in a browser for executing OncoBird. If you want to build your own container, run `./environment/make_docker.sh`

Data input format
--
The molecular and the clinical data for running OncoBird can be jointly in one data frame, or supplied as clinical data and mutational data separately. Thereby, the following format is expected.
```
data/data_mutations.csv
```
| sample | \<MUT1\>_AMP | \<MUT2\>_SV | \<MUT3\>_DEL |
| ------ | ----------- | ----------- | ----------- |
| \<X1\> | 1 | 1 | 1 |
| \<X2\> | 1 | 0 | 0 |
| \<X3\> | 0 | 1 | 0 |
| \<X4\> | 0 | 0 | 0 |
| ... |  ... | ...  | ...  |

```
data/data_clinical.csv
```
| sample | treatment | \<subtype1\> | \<subtype2\> | \<OS/PFS\>.event | \<OS/PFS\>.months | ORR |
| ------ | ----------- | ----------- |----------- |----------- |----------- |----------- |
| \<X1\> | t1 | C1 | left | 0| 11 | 0 |
| \<X2\> | t1 | C2 | right | 0| 22| 1|
| \<X3\> | t2 | C2 | left | 0| 33| 1|
| \<X4\> | t2 | C2 | left | 1| 44| 1|
| ...| | ...  | ...  | ...  | ... | ...| 


Reference
--
For citing OncoBird, please refer to the original manuscript [DOI pending].
