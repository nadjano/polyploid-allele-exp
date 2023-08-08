

### Installing Nextflow :

```
wget -qO- https://get.nextflow.io | bash

```


#### Running nextflow command

```
./nextflow run file.nf
```


### Information

#### WHy BWA ?
Our reads are long (250 bp) so we will use BWA-MEM (Li 2013) to align them against the reference genome as it has good mapping performance for longer reads (100bp and up)

