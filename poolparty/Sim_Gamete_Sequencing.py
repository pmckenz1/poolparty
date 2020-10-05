import os
import numpy as np
import h5py
from .utils import *
from tqdm import tqdm_notebook


class Sim_Gamete_Sequencing():
    def __init__(
        self,
        directory,
        pdf,
        num_gams,  # total number of gametes
        gpa=None,  # gams per aliquot (this x nali should equal num_gams)
        nali=None,  # num aliquots (this x gpa should equal num_gams)
        ncutsites=None,  # num cut sites per genome
        num_reads=None,  # num reads to generate
        directory_exists=True,
    ):
        self.num_gams = num_gams
        self.gpa = gpa
        self.nali = nali
        self.ncutsites = ncutsites
        self.num_reads = num_reads
        self.directory = directory
        self.directory_exists = True

        self.pdf = pdf

    def sim_gametes_only(self, aliquots=False):

        if not self.directory_exists:
            os.mkdir(self.directory)

        self.gamspath = os.path.join(self.directory, 'gams.hdf5')

        gams = h5py.File(self.gamspath,'w')

        gams.create_dataset("start_haplo",
                            data=np.random.binomial(1,.5,size=self.num_gams),
                            dtype=np.int8)

        gams.create_dataset("num_crossovers",
                            data=1+np.random.binomial(1,.2,size=self.num_gams),
                            dtype=np.int8)

        np.sum(gams["num_crossovers"])

        gams.create_dataset("crossovers",
                            data = np.vstack([np.hstack([np.repeat(i,gams['num_crossovers'][i]) for i in range(self.num_gams)]),
                                              self.pdf.rvs(size=np.sum(gams["num_crossovers"]))]).T)

        gams.close()
        
        if aliquots:
            self.aliquotspath = os.path.join(self.directory, 'aliquots.hdf5')
            gams = h5py.File(self.gamspath,'r')

            samp_arr = np.random.choice(range(gams['num_crossovers'].shape[0]),
                             size=self.gpa*self.nali,
                             replace=False).reshape((self.nali,self.gpa))

            with h5py.File(self.aliquotspath,"w") as f:
                f.create_dataset("aliquots",
                                data=samp_arr,
                                dtype=np.int64)

            gams.close()
        
    def sim_gametes_and_sequencing(self,
                                   evenly_spaced_loci=False
                                   ):
        print("Simulating gametes...")
        self.sim_gametes_only(aliquots=True)
        
        samp_arr_file = h5py.File(self.aliquotspath, "r")

        samp_arr = samp_arr_file['aliquots']
        
        print("Sequencing gametes...")
        self.seqspath = os.path.join(self.directory, 'seqs.hdf5')
        with h5py.File(self.seqspath,'w') as f:
            f.create_dataset("seqs",
                            data=do_sequencing(samp_arr, 
                                               self.ncutsites, 
                                               self.num_reads),
                            dtype=np.int64)

        seqh5 = h5py.File(self.seqspath, 'r')
        
        print("Demultiplexing...")
        self.seqs_reshapepath = os.path.join(self.directory, 'seqs_reshape.hdf5')
        seqs_reshape = h5py.File(self.seqs_reshapepath,'w')
        for alinum in tqdm_notebook(range(self.nali)):
            tmpgroup = seqs_reshape.create_group(str(alinum))
            for locusnum in range(self.ncutsites):
                alisamp = np.take(seqh5['seqs'],np.where(seqh5['seqs'][:,0] == alinum),axis=0)[0]
                gams_seqed = np.take(alisamp,np.where(alisamp[:,2] == locusnum),axis=0)[0][:,1]

                tmpgroup.create_dataset(str(locusnum),
                                  dtype=np.int64,
                                  data=gams_seqed)

        seqs_reshape.close()

        # now start the haplotyping, which is super slow....
        
        print('Assigning haplotypes to loci -- this might take a while...')
        
        gams = h5py.File(self.gamspath,'r')

        crossarr = gams['crossovers']
        haplostart = gams['start_haplo']

        seqs_reshape = h5py.File(self.seqs_reshapepath,'r')
        
        self.haplotypespath = os.path.join(self.directory, 'haplotypes.hdf5')
        haplotypes = h5py.File(self.haplotypespath,'w')

        if not evenly_spaced_loci:
            loci = np.sort(np.random.uniform(0,1,self.ncutsites))
        else:
            loci = np.linspace(0,1,self.ncutsites)
        haplotypes.create_dataset("loci_locs",
                                  data=loci
                                  )

        for ali_idx_ in tqdm_notebook(range(self.nali)):
            tmpgroup = haplotypes.create_group(str(ali_idx_))
            for locus_idx_ in range(self.ncutsites):

                readidxs = seqs_reshape[str(ali_idx_)][str(locus_idx_)]
                readloc = loci[locus_idx_]
                haplos = []
                for readidx in readidxs:
                    crosses = np.take(crossarr[:,1],np.where(crossarr[:,0] == readidx))[0]
                    #crosses = crossarr[crossarr[:,0] == readidx,1]
                    starter = np.take(haplostart,
                                      readidx)

                    passes = np.sum(readloc > crosses)
                    #read_haplo = starter
                    #for i in range(passes):
                    #    read_haplo = 1-read_haplo

                    haplos.append(convert_haplo(passes,starter))

                tmpgroup.create_dataset(str(locus_idx_),
                                  dtype=np.int64,
                                  data=haplos)
            #print(ali_idx_)

        haplotypes.close()


def do_sequencing(arr,
                  ncutsites,
                  num_reads,
                 ):
    seqarr = np.zeros((num_reads,3),dtype=np.int64)

    for rowidx in tqdm_notebook(range(num_reads)):
        seqarr[rowidx] = get_read(arr,ncutsites)
    return(seqarr)