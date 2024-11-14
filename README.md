## long-read ASE


Prequesites create allelefinder env


# Install AlleleFinder
mamba create -n allelefinder gmap blast
conda activate allelefinder
cd /scratch/nadjafn
git clone https://github.com/sc-zhang/AlleleFinder.git
chmod +x AlleleFinder/allelefinder.py
# Optional
echo 'export PATH=/scratch/nadjafn/AlleleFinder:$PATH' >> ~/.bash_profile
echo 'export PATH=/scratch/nadjafn/MCScanX/MCScanX:$PATH' >> ~/.bash_profile
source ~/.bash_profile

# nf-potato-ase
