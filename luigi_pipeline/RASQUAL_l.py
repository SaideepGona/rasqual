'''
Author : Saideep Gona

This PYTHON3 script is intended for running an end-to-end 
RASQUAL pipeline with any number of data modalities.

Based off the code from: https://github.com/natsuhiko/rasqual
With modification to improve scalability and flexibility on 
the cluster.

'''

import sys, os, argparse
from glob import glob
from datetime import datetime
from os.path import join

import linecache

import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi.util import requires
from bioluigi.scheduled_external_program import ScheduledExternalProgramTask

import threading
import subprocess
import time

# Boilerplate
cwd = os.getcwd()           # Store current working directory
print("Current Working Directory: ", cwd)
dtime = str(datetime.now())
print("Time of Running ", dtime)
print("Python Interpreter: ", sys.executable)

class CustomWrapperTask(luigi.WrapperTask):
    CACHED_REQUIRES = []

    def cached_requires(self):
        return self.CACHED_REQUIRES or self.requires()

    def complete(self):
        return all(r.complete() for r in self.cached_requires())

    def deps(self):
        return self.cached_requires()

#####################################################################################
# CONFIG AND PREP INPUT FILES
#####################################################################################

class RASQUAL(luigi.Config):

    GENOME = luigi.Parameter()
    SPLIT_GENOME_DIR = luigi.Parameter()
    CHROM_INFO = luigi.Parameter()
    STAR_INDEX = luigi.Parameter()
    ANNOTATION = luigi.Parameter()

    PHASING_TYPE = luigi.Parameter()
    PHASED_DATA_DIR = luigi.Parameter()

    PERSONAL_VCF = luigi.Parameter()

    SAMPLE_ORDER = luigi.Parameter()

    SINGLE_END = luigi.Parameter()

    ASSAY = luigi.Parameter()

    USER = luigi.Parameter()

    WINDOW_SIZE = luigi.IntParameter()

    OUTPUT_DIR = luigi.Parameter()

    sbatch_templates = luigi.Parameter()

    cluster_max_jobs = luigi.IntParameter()

    samtools_bin = luigi.Parameter()
    star_bin = luigi.Parameter()
    gatk_bin = luigi.Parameter()
    tabix_bin = luigi.Parameter()
    filter_wasp_path = luigi.Parameter()
    rasqual_dir = luigi.Parameter()
    create_rasqual_tests_path = luigi.Parameter()
    stringtie_path = luigi.Parameter()
    prepDE_path = luigi.Parameter()
    rasqual_tools_path = luigi.Parameter()

cfg = RASQUAL()

os.system("mkdir " + cfg.OUTPUT_DIR + "/logs")

class ProducePhasedVCF(luigi.Task):
    """
    Produce an annotated VCF file
    """
    def output(self):
        return luigi.LocalTarget(cfg.PERSONAL_VCF)

class ProducePhasedVCFGZ(luigi.Task):
    """
    Produce an annotated VCF file
    """
    def output(self):
        return luigi.LocalTarget(cfg.PERSONAL_VCF+".gz")

class ProducePhasedVCFGZ_Tabix(luigi.Task):
    """
    Produce an annotated VCF file
    """
    def output(self):
        return luigi.LocalTarget(cfg.PERSONAL_VCF+".gz.tbi")

class ProduceStarIndex(luigi.Task):

    def output(self):
        return luigi.LocalTarget(join(cfg.STAR_INDEX,"SA"))

class ProduceGenome(luigi.Task):
    """
    Produce a reference genome sequence.
    """
    def output(self):
        return luigi.LocalTarget(cfg.GENOME)

class ProduceAnnotations(luigi.Task):
    """
    Produce a reference annotations
    """
    def output(self):
        return luigi.LocalTarget(cfg.ANNOTATION)

class ProduceSampleFastqs(luigi.Task):
    """
    Produce the FASTQs that relate to a sample.
    """
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()
    def output(self):
        return [luigi.LocalTarget(f) for f in sorted(glob(join(cfg.OUTPUT_DIR, 'data', self.experiment_id, self.sample_id, '*.fastq.gz')))]

@requires(ProducePhasedVCF)
class InputOrder(luigi.Task):
    experiment_id = luigi.Parameter()
    def run(self):
        with open(self.input().path, "r") as vcf:
            for line in vcf:
                if line.startswith("#CHROM"):
                    p_line = line.split("\t")
                    print(p_line)
                    line_len = len(p_line)
                    if "FORMAT" in line:
                        num_samples = line_len - 9
                    else:
                        num_samples = line_len - 8
        with open(cfg.SAMPLE_ORDER, "r") as samples:
            sr = samples.readlines()
            print(len(sr))
            print(num_samples)
            if len(sr) == num_samples:
                os.system("touch " + self.output().path)
            else:
                raise Exception("Number of samples in sample ordering file and vcf headers don't match")

    def output(self):
        return luigi.LocalTarget(cfg.SAMPLE_ORDER +".confirmed")

#####################################################################################
# WASP ALIGNMENT
#####################################################################################

# @requires(ProduceStarIndex, ProducePhasedVCF, ProduceSampleFastqs)
# class AlignSampleWASP(ScheduledExternalProgramTask):
#     """
#     Align a sample on a reference.
#     """
#     scheduler = 'slurm'
#     cpus = 8
#     memory = 50

#     def program_args(self):

#         args = [       
#             cfg.star_bin, 
#             "--outFileNamePrefix", os.path.dirname(self.output().path) + '/',
#             "--runThreadN", self.cpus,
#             "--genomeDir", cfg.STAR_INDEX,
#             "--outSAMstrandField", "intronMotif",
#             "--outSAMtype", "BAM SortedByCoordinate",
#             "--outFilterIntronMotifs", "RemoveNoncanonical",
#             "--varVCFfile", self.input()[1].path,
#             "--waspOutputMode", "SAMtag",
#             "--outSAMattributes", "vA vG",
#             "--readFilesCommand", "zcat"
#         ]

#         args.append('--readFilesIn')
#         args.extend(fastq.path for fastq in self.input()[2])

#         return args

#     def run(self):
#         self.output().makedirs()
#         return super(AlignSampleWASP, self).run()

#     def output(self):
#         return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'aligned', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.bam'))

# @requires(AlignSampleWASP)
# class IndexAlignment(ScheduledExternalProgramTask):
#     """
#     Perform an alignment (2-step)
#     """
#     scheduler = 'slurm'
#     resources = {'cpu': 8, 'mem': 32}

#     def program_args(self):
#         args = [cfg.samtools_bin, "index",
#                 self.input().path,
#                 self.output().path]

#         return args

#     def run(self):
#         self.output().makedirs()
#         return super(IndexAlignment, self).run()

#     def output(self):
#         return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'aligned-step2', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.bai'))


# @requires(IndexAlignment)
# class FilterWASP(ScheduledExternalProgramTask):
#     """
#     Perform an alignment (2-step)
#     """
#     scheduler = 'slurm'
#     resources = {'cpu': 1, 'mem': 2}

#     def program_args(self):
#         args = [cfg.filter_wasp_path,
#                 self.input().path,
#                 self.output().path]

#         return args

#     def run(self):
#         self.output().makedirs()
#         return super(FilterWASP, self).run()

#     def output(self):
#         return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'filter_wasp', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.filtered.bam'))

# @requires(FilterWASP)
# class IndexAlignmentFiltered(ScheduledExternalProgramTask):
#     """
#     Perform an alignment (2-step)
#     """
#     scheduler = 'slurm'
#     resources = {'cpu': 8, 'mem': 32}

#     def program_args(self):
#         args = [cfg.samtools_bin, "index",
#                 self.input().path,
#                 self.output().path]

#         return args

#     def run(self):
#         self.output().makedirs()
#         return super(IndexAlignmentFiltered, self).run()

#     def output(self):
#         return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'filter_wasp', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.filtered.bai'))

# class WASP_alignment(luigi.Task):

#     experiment_id = luigi.Parameter()

#     # resources = {'cpu': 8, 'mem': 32}

#     def run(self):
#         samples = [sample_id for sample_id in glob(join(cfg.OUTPUT_DIR, 'data', self.experiment_id, '*'))]
        
#         # chunk = 3
#         # for x in [1,2]:
#         #     procs = []
#         for sample in samples[(x-1)*chunk:x*chunk]:
#             print(sample)
#             command = [
#                 "bash",
#                 "luigi-wrapper",
#                 "IndexAlignmentFiltered",
#                 "--experiment-id", self.experiment_id,
#                 "--sample-id", sample
#             ]


# class CheckPoint1(CustomWrapperTask):
#     experiment_id = luigi.Parameter()
#     def requires(self):
#         reqs = [IndexAlignmentFiltered(self.experiment_id, os.path.basename(sample_id))
#                 for sample_id in glob(join(cfg.OUTPUT_DIR, 'data', self.experiment_id, '*'))]
#         for req in reqs:
#             self.CACHED_REQUIRES.append(req)
#         return self.CACHED_REQUIRES

#####################################################################################
# WASP ALIGNMENT (STANDARD)
#####################################################################################

@requires(ProduceStarIndex, ProducePhasedVCF, ProduceSampleFastqs)
class AlignSampleWASP(ExternalProgramTask):
    """
    Align a sample on a reference.
    """

    def program_args(self):

        args = [       
            cfg.star_bin, 
            "--outFileNamePrefix", os.path.dirname(self.output().path) + '/',
            "--runThreadN", "12",
            "--genomeDir", cfg.STAR_INDEX,
            "--outSAMstrandField", "intronMotif",
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outFilterIntronMotifs", "RemoveNoncanonical",
            "--varVCFfile", self.input()[1].path,
            "--waspOutputMode", "SAMtag",
            "--outSAMattributes", "vA", "vG",
            "--readFilesCommand", "zcat"
        ]

        args.append('--readFilesIn')
        args.extend(fastq.path for fastq in self.input()[2])

        return args

    def run(self):
        self.output().makedirs()
        return super(AlignSampleWASP, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'aligned', self.experiment_id, self.sample_id, 'Log.final.out'))

@requires(AlignSampleWASP)
class IndexAlignment(ExternalProgramTask):
    """
    Perform an alignment (2-step)
    """

    def program_args(self):
        args = [cfg.samtools_bin, "index",
                self.input().path,
                self.output().path]

        return args

    def run(self):
        self.output().makedirs()
        return super(IndexAlignment, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'aligned', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.bai'))


@requires(IndexAlignment)
class FilterWASP(ExternalProgramTask):
    """
    Perform an alignment (2-step)
    """

    def program_args(self):
        args = [cfg.filter_wasp_path,
                self.input().path,
                self.output().path]

        return args

    def run(self):
        self.output().makedirs()
        return super(FilterWASP, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'filter_wasp', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.filtered.bam'))

@requires(FilterWASP)
class IndexAlignmentFiltered(ExternalProgramTask):
    """
    Perform an alignment (2-step)
    """

    def program_args(self):
        args = [cfg.samtools_bin, "index",
                self.input().path,
                self.output().path]

        return args

    def run(self):
        self.output().makedirs()
        return super(IndexAlignmentFiltered, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'filter_wasp', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.filtered.bai'))

@requires(InputOrder)
class WASP_alignment(luigi.Task):

    experiment_id = luigi.Parameter()

    def run(self):
        samples = []
        with open(cfg.SAMPLE_ORDER, "r") as inp:
            for line in inp:
                samples.append(line.rstrip("\n"))
        
        # samples = [sample_id for sample_id in glob(join(cfg.OUTPUT_DIR, 'data', self.experiment_id, '*'))]

        for sample in samples:
            command = [
                "bash",
                "luigi-wrapper",
                "IndexAlignmentFiltered",
                "--experiment-id", self.experiment_id,
                "--sample-id", sample
            ]
            run_file = join(cfg.OUTPUT_DIR, "logs", "WASP_align_" + self.experiment_id + "_" + sample + ".sbatch")
            os.system("cp " + join(cfg.sbatch_templates,"WASP_align.sbatch") + " " + run_file)
            with open(run_file, "a") as runsb:
                runsb.write("#SBATCH --job-name=WASP_align_" + sample + "\n")
                runsb.write("#SBATCH -o "+join(cfg.OUTPUT_DIR, "logs", "WASP_align_" + self.experiment_id + "_" + sample + "%J.o")+"\n")
                runsb.write("#SBATCH -e "+join(cfg.OUTPUT_DIR, "logs", "WASP_align_" + self.experiment_id + "_" + sample + "%J.e")+"\n")
                runsb.write("\n" + " ".join(command) + "\n")
            os.system("sbatch "  + run_file)

    def output(self):

        samples = []
        with open(cfg.SAMPLE_ORDER, "r") as inp:
            for line in inp:
                samples.append(line.rstrip("\n"))

        return [luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'filter_wasp', self.experiment_id, sample, 'Aligned.sortedByCoord.out.filtered.bai')) for sample in samples]


@requires(InputOrder)
class Check_WASP_alignment(luigi.WrapperTask):
    experiment_id = luigi.Parameter()
    def requires(self):
        samples = []
        with open(cfg.SAMPLE_ORDER, "r") as inp:
            for line in inp:
                samples.append(line.rstrip("\n"))
        return [IndexAlignmentFiltered(self.experiment_id, sample) for sample in samples]
        



#####################################################################################
# STRINGTIE ASSEMBLY
#####################################################################################


@requires(ProduceAnnotations, IndexAlignmentFiltered)
class StringtieAssembly(ScheduledExternalProgramTask):
    """
    
    """
    scheduler = 'slurm'
    resources = {'cpu': 8, 'mem': 16}

    def program_args(self):
        args = [
            cfg.stringtie_path,
            self.input()[1].path,
            "-o", self.output().path,
            "-p", "8",
            "-e",
            "-G", self.input()[0].path
            ]

        return args

    def run(self):
        self.output().makedirs()
        return super(StringtieAssembly, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'stringtie', self.experiment_id, self.sample_id, 'assembly.gtf'))

class StringtieMergeList(luigi.Task):
    
    def requires(self):
        require = [InputOrder(self.experiment_id, self.sample_id)]
        assemblies = [StringtieAssembly(self.experiment_id, os.path.basename(sample_id))
                for sample_id in glob(join(cfg.OUTPUT_DIR, 'data', self.experiment_id, '*'))]
        return require + assemblies
        
    def run(self):
        self.output().makedirs()
        with open(self.input().path, "r") as inp:
            samples = inp.readlines()

        with open(self.output().path, "w") as o:
            for sample in samples:
                o.write(join(cfg.OUTPUT_DIR, 'stringtie', self.experiment_id, sample, 'assembly.gtf')+"\n")


    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'stringtie', self.experiment_id, 'merge_list.txt'))

@requires(ProduceAnnotations, StringtieMergeList)
class StringtieMerge(ScheduledExternalProgramTask):

    scheduler = 'slurm'
    resources = {'cpu': 8, 'mem': 16}


    def program_args(self):

        args = [
            cfg.stringtie_path,
            "--merge",
            "-p", "8",
            "-G", self.input()[0].path,
            "-o", self.output().path,
            "-e", 
            self.input()[1].path
        ]

        return args

    def run(self):
        self.output().makedirs()
        return super(StringtieMerge, self).run()      

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'stringtie', self.experiment_id, 'merged.gtf'))

@requires(StringtieMerge)
class MSTRGPrep(ScheduledExternalProgramTask):

    scheduler = 'slurm'
    resources = {'cpu': 1, 'mem': 2}

    def program_args(self):

        args = [
            cfg.mstrg_prep_path,
            self.input().path,
            ">",
            self.output().path
        ]

        return args

    def run(self):
        self.output().makedirs()
        return super(StringtieMerge, self).run()      

    def output(self):
        return luigi.LocalTarget(join(self.input().path.rstrip(".gtf"), "_mstrgprep.gtf"))


@requires(MSTRGPrep, IndexAlignmentFiltered)
class FinalAbundance(ScheduledExternalProgramTask):

    scheduler = 'slurm'
    resources = {'cpu': 8, 'mem': 16}

    def program_args(self):

        args = [
            cfg.stringtie_path,
            self.input()[1].path,
            "-e",
            "-p", "8",
            "-G", self.input()[0].path,
            "-o", self.output().path
        ]

        return args

    def run(self):

        self.output().makedirs()
        return super(FinalAbundance, self).run()      

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'stringtie_final_output', self.experiment_id, self.sample_id, 'final.gtf'))

class CreateAbundanceFile(luigi.Task):

    def requires(self):

        require = [InputOrder(self.experiment_id, self.sample_id)]

        with open(self.input().path, "r") as inp:
            samples = inp.readlines()

        with open(self.output().path, "w") as o:
            for sample in samples:
                require += FinalAbundance(self.experiment_id, os.path.basename(sample))

        return require

    def run(self):
        self.output().makedirs()
        with open(self.input().path, "r") as inp:
            samples = inp.readlines()

        with open(self.output().path, "w") as o:
            for sample in samples:
                o.write(join(cfg.OUTPUT_DIR, 'stringtie_final_output', self.experiment_id, sample, 'final.gtf')+"\n")

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'final_counts', self.experiment_id, 'stringtie_counts.txt'))

@requires(CreateAbundanceFile)
class GenerateCounts(ScheduledExternalProgramTask):

    scheduler = 'slurm'
    resources = {'cpu': 1, 'mem': 2}

    def program_args(self):

        args = [
            cfg.prepDE_path,
        "-i", self.input().path,
        "-g", join(cfg.OUTPUT_DIR, 'final_counts', self.experiment_id, 'gene_counts.csv'),
        "-t", join(cfg.OUTPUT_DIR, 'final_counts', self.experiment_id, 'transcript_counts.csv')
        ]

        return args

    def run(self):

        self.output().makedirs()
        return super(GenerateCounts, self).run()      

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'final_counts', self.experiment_id, 'gene_counts.csv'))

@requires(GenerateCounts)
class PrepExpression(ScheduledExternalProgramTask):

    scheduler = 'slurm'
    resources = {'cpu': 1, 'mem': 2}

    def program_args(self):

        args = [
            cfg.prep_expressionset_path,
            self.input().path,
            self.output().path
        ]
        return args

    def run(self):
        self.output().makedirs()
        return super(PrepExpression, self).run()      

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, 'final_counts', self.experiment_id, 'assay_data.csv'))



#####################################################################################
# RASQUALTOOLS
#####################################################################################


@requires(GenerateCounts)
class RasqualTools(ScheduledExternalProgramTask):

    scheduler = 'slurm'
    resources = {'cpu': 2, 'mem': 10}

    def program_args(self):

        args = [
            "Rscript",
            cfg.rasqual_tools_path + "/rasqualTools.R", 
            self.input().path,
            os.path.dirname(self.output()[0].path)
        ]
        return args

    def run(self):
        self.output().makedirs()
        return super(RasqualTools, self).run()      

    def output(self):
        return [luigi.LocalTarget(join(cfg.OUTPUT_DIR,"binary_read_counts",self.experiment_id,"Y.bin")),
        luigi.LocalTarget(join(cfg.OUTPUT_DIR,"binary_read_counts",self.experiment_id,"K.bin"))]

class Raw_Read_Counts(luigi.Task):

    experiment_id = luigi.Parameter()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR,"final_counts",self.experiment_id,"assay_data.csv"))


class Read_Counts(luigi.Task):

    experiment_id = luigi.Parameter()

    def output(self):
        return [luigi.LocalTarget(join(cfg.OUTPUT_DIR,"binary_read_counts",self.experiment_id,"Y.bin")),
        luigi.LocalTarget(join(cfg.OUTPUT_DIR,"binary_read_counts",self.experiment_id,"K.bin"))]

#####################################################################################
# CREATE VCF WITH ALLELE SPECIFIC COUNTS
#####################################################################################


class Create_Bam_List(luigi.Task):

    experiment_id = luigi.Parameter()

    def requires(self):
        return [AlignSampleWASP(self.experiment_id, os.path.basename(sample_id))
                for sample_id in glob(join(cfg.OUTPUT_DIR, 'data', self.experiment_id, '*'))]

    def run(self):
        self.output().makedirs()
        print(self.input())
        with open(self.output().path, "w") as out:
            sort_in = [x.path for x in self.input()]
            sort_in.sort()
            for bam in sort_in:
                out.write(bam + "\n")

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "allele_specific", self.experiment_id, "bam_list.txt"))

@requires(ProducePhasedVCFGZ, ProducePhasedVCFGZ_Tabix, Create_Bam_List)
class ASVCF(ExternalProgramTask):
    """
    
    """
    # scheduler = 'slurm'
    # resources = {'cpu': 8, 'mem': 32}

    def program_args(self):
        os.system("export RASQUALDIR=" + cfg.rasqual_dir)
        if cfg.SINGLE_END == "True":
            end = "single_end"
        else:
            end = "paired_end"
        args = [
            "bash",
            join(cfg.rasqual_dir,"src","ASVCF","createASVCF.sh"),
            end,
            self.input()[2].path,
            self.input()[0].path,
            self.output().path,
            cfg.ASSAY
            ]
        return args

    def run(self):
        self.output().makedirs()
        return super(ASVCF, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "allele_specific", self.experiment_id, "ASVCF.vcf.gz"))

@requires(ASVCF)
class ASVCF_Tabix(ExternalProgramTask):

    def program_args(self):

        args = [
            cfg.tabix_bin,
            "-p", "vcf",
            self.input().path
            ]
        return args

    def run(self):
        self.output().makedirs()
        return super(ASVCF_Tabix, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "allele_specific", self.experiment_id, "ASVCF.vcf.gz.tbi"))


#####################################################################################
# GENERATE RASQUAL INPUT PARAMS
#####################################################################################

@requires(ProduceAnnotations)
class Create_Windowed_GTF(luigi.Task):
    """
    Creates a new gtf file with gene start/end expanded to the window size
    """
    
    experiment_id = luigi.Parameter()

    def run(self):
        print(self.input().path)
        
        with open(self.input().path, "r") as inp:
            script = ""    
            for line in inp:
                if line.startswith("#"):
                    continue
                p_line = line.split("\t")
                start = p_line[3]
                new_start = int(p_line[3]) - cfg.WINDOW_SIZE
                if new_start < 0:
                    new_start = 0
                new_start = str(new_start)
                end = p_line[4]
                new_end = str(int(p_line[4]) + cfg.WINDOW_SIZE)
                new_line = p_line[:]
                new_line[3] = new_start
                new_line[4] = new_end
                with open(self.output().path, "a") as out:
                    out.write("\t".join(new_line) + "\n")


    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "generate_rasqual_tests", self.experiment_id, "windowed.gtf"))

@requires(Create_Windowed_GTF, ProducePhasedVCF)
class Intersect_GTF_VCF(luigi.Task):
    """
    Creates a new gtf file with gene start/end expanded to the window size
    """

    def run(self):
        command = [
            "bedtools intersect",
            "-a", self.input()[0].path,
            "-b", self.input()[1].path,
            "-wa", "-wb",
            ">",
            self.output().path
        ]
        os.system(" ".join(command))

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "generate_rasqual_tests", self.experiment_id, "intersect_gtf_vcf.gtf"))

@requires(Intersect_GTF_VCF, Raw_Read_Counts)
class Create_RASQUAL_Tests(ScheduledExternalProgramTask):
    """
    Creates semi-populated RASQUAL commands ready for execution
    """
    scheduler = 'slurm'
    cpus = 8
    memory = 30



    def program_args(self):

        args = [
            "python",
            cfg.create_rasqual_tests_path,
            self.input()[0].path,
            self.output().path,
            self.input()[1].path
            ]
        return args

    def run(self):
        self.output().makedirs()
        return super(Create_RASQUAL_Tests, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "generate_rasqual_tests", self.experiment_id, "RASQUAL_tests.txt"))



#####################################################################################
# RUN RASQUAL
#####################################################################################

# @requires(ASVCF_Tabix, Read_Counts, Create_RASQUAL_Tests)
# class Run_RASQUAL(ScheduledExternalProgramTask):
#     """
#     Creates semi-populated RASQUAL commands ready for execution
#     """
#     scheduler = 'slurm'
#     cpus = 1
#     memory = 1

#     experiment_id = luigi.Parameter()

#     feature_id = luigi.Parameter()


#     def program_args(self):

#         script = linecache.getline(self.input()[2].path, int(self.feature_id))
            
#         script = script.replace("PLACEHOLDER_TABIX", cfg.tabix_bin)
#         script = script.replace("PLACEHOLDER_VCF", self.input()[0].path.rstrip(".tbi"))
#         script = script.replace("PLACEHOLDER_RASQUAL", join(cfg.rasqual_dir,"bin","rasqual"))
#         script = script.replace("PLACEHOLDER_VCF", self.input()[0].path)
#         script = script.replace("PLACEHOLDER_YBIN", self.input()[1][0].path)
#         script = script.replace("PLACEHOLDER_KBIN", self.input()[1][1].path)
#         script = script.replace("PLACEHOLDER_NUMTHREADS", "10")

#         args = script.split(" ")
#         args.append(">")
#         args.append(self.output().path)

#         return args

#     def run(self):
#         self.output().makedirs()
#         return super(Run_RASQUAL, self).run()

#     def output(self):
#         return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "rasqual_out", self.experiment_id,str(self.feature_id) + ".out"))

@requires(ASVCF_Tabix, PrepExpression, Create_RASQUAL_Tests)
class Run_RASQUAL(luigi.Task):
    """
    Runs a RASQUAL task on a single feature
    """
    # scheduler = 'slurm'
    # cpus = 1
    # memory = 1

    experiment_id = luigi.Parameter()

    feature_id = luigi.Parameter()


    def run(self):
        self.output().makedirs()

        script = linecache.getline(self.input()[2].path, int(self.feature_id)).rstrip("\n")
            
        script = script.replace("PLACEHOLDER_TABIX", cfg.tabix_bin)
        script = script.replace("PLACEHOLDER_VCF", self.input()[0].path.rstrip(".tbi"))
        script = script.replace("PLACEHOLDER_RASQUAL", join(cfg.rasqual_dir,"bin","rasqual"))
        script = script.replace("PLACEHOLDER_VCF", self.input()[0].path)
        script = script.replace("PLACEHOLDER_YBIN", self.input()[1][0].path)
        script = script.replace("PLACEHOLDER_KBIN", self.input()[1][1].path)
        script = script.replace("PLACEHOLDER_NUMTHREADS", "10")

        args = script.split(" ")
        args.append(">"),
        args.append(self.output().path)

        subprocess.run(" ".join(args), shell=True)

        run_file = join(cfg.OUTPUT_DIR, "logs",self.experiment_id, "RASQUAL_" + str(self.feature_id) + ".run")
        with open (run_file, "w") as r:
            r.write(" ".join(args))

        # args = script.split("|")
        # all_args = [x.split(" ") for x in args]
        # r1 = subprocess.run(" ".join(all_args[0]),text=True, 
        # shell=True, stdout=subprocess.PIPE)
        # r2 = subprocess.run(" ".join(all_args[1]), shell=True,
        #     stdin = r1.stdout, text=True, stdout=self.output().open(mode="w"))
        return super(Run_RASQUAL, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "rasqual_out", self.experiment_id,str(self.feature_id) + ".out"))



@requires(Create_RASQUAL_Tests)
class Run_All_RASQUAL(luigi.Task):

    experiment_id = luigi.Parameter()

    def run(self):
        
        with open(self.input().path, "r") as tests:
            all_tests = tests.readlines()

        submits = 1
        for i,test in enumerate(all_tests):
            if submits > cfg.cluster_max_jobs:
                break
            if os.path.exists(join(cfg.OUTPUT_DIR, "rasqual_out", self.experiment_id,str(i + 1) + ".out")):
                continue
            submits += 1
            command = [
                "bash",
                "luigi-wrapper",
                "Run_RASQUAL",
                "--experiment-id", self.experiment_id,
                "--feature-id", str(i+1)
            ]
            os.system("mkdir -p " + join(cfg.OUTPUT_DIR, "logs",self.experiment_id))
            run_file = join(cfg.OUTPUT_DIR, "logs",self.experiment_id, "RASQUAL_" + str(i) + ".sbatch")
            
            os.system("cp " + join(cfg.sbatch_templates,"RASQUAL.sbatch") + " " + run_file)
            with open(run_file, "a") as runsb:
                runsb.write("#SBATCH --job-name=RASQUAL_" + self.experiment_id + "_" + str(i) + "\n")
                runsb.write("#SBATCH -o "+join(cfg.OUTPUT_DIR, "logs", self.experiment_id, "RASQUAL_" + str(i) + "_" + "%J.o")+"\n")
                runsb.write("#SBATCH -e "+join(cfg.OUTPUT_DIR, "logs", self.experiment_id, "RASQUAL_" + str(i) + "_" + "%J.e")+"\n")
                runsb.write("\n" + " ".join(command) + "\n")
            os.system("sbatch " + run_file)

class Check_All_RASQUAL(luigi.Task):

    experiment_id = luigi.Parameter()

    def run(self):

        tests = Create_RASQUAL_Tests(self.experiment_id)
        yield tests

        with open(tests.output().path, "r") as t:
            all_tests = t.readlines()

        for i,test in enumerate(all_tests):
            yield Run_RASQUAL(self.experiment_id,str(i))

#####################################################################################
# RUN RASQUAL PERMUTATIONS
#####################################################################################

@requires(ASVCF_Tabix, PrepExpression, Create_RASQUAL_Tests)
class Run_RASQUAL_permutation(luigi.Task):
    """
    Runs a RASQUAL task on a single feature
    """
    # scheduler = 'slurm'
    # cpus = 1
    # memory = 1

    experiment_id = luigi.Parameter()

    feature_id = luigi.Parameter()


    def run(self):
        self.output().makedirs()

        script = linecache.getline(self.input()[2].path, int(self.feature_id)).rstrip("\n")
            
        script = script.replace("PLACEHOLDER_TABIX", cfg.tabix_bin)
        script = script.replace("PLACEHOLDER_VCF", self.input()[0].path.rstrip(".tbi"))
        script = script.replace("PLACEHOLDER_RASQUAL", join(cfg.rasqual_dir,"bin","rasqual"))
        script = script.replace("PLACEHOLDER_VCF", self.input()[0].path)
        script = script.replace("PLACEHOLDER_YBIN", self.input()[1][0].path)
        script = script.replace("PLACEHOLDER_KBIN", self.input()[1][1].path)
        script = script.replace("PLACEHOLDER_NUMTHREADS", "10")

        args = script.split(" ")
        args.append("-r")
        args.append(">"),
        args.append(self.output().path)

        subprocess.run(" ".join(args), shell=True)

        run_file = join(cfg.OUTPUT_DIR, "logs",self.experiment_id, "RASQUAL_" + str(self.feature_id) + ".run")
        with open (run_file, "w") as r:
            r.write(" ".join(args))


        return super(Run_RASQUAL, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.OUTPUT_DIR, "rasqual_out", self.experiment_id,str(self.feature_id) + ".out"))



@requires(Create_RASQUAL_Tests)
class Run_All_RASQUAL_permutations(luigi.Task):

    experiment_id = luigi.Parameter()

    def run(self):
        
        with open(self.input().path, "r") as tests:
            all_tests = tests.readlines()

        submits = 1
        for i,test in enumerate(all_tests):
            if submits > cfg.cluster_max_jobs:
                break
            if os.path.exists(join(cfg.OUTPUT_DIR, "rasqual_out", self.experiment_id,str(i + 1) + ".out")):
                continue
            submits += 1
            command = [
                "bash",
                "luigi-wrapper",
                "Run_RASQUAL_permutation",
                "--experiment-id", self.experiment_id,
                "--feature-id", str(i+1)
            ]
            os.system("mkdir -p " + join(cfg.OUTPUT_DIR, "logs",self.experiment_id))
            run_file = join(cfg.OUTPUT_DIR, "logs",self.experiment_id, "RASQUAL_PERM_" + str(i) + ".sbatch")
            
            os.system("cp " + join(cfg.sbatch_templates,"RASQUAL.sbatch") + " " + run_file)
            with open(run_file, "a") as runsb:
                runsb.write("#SBATCH --job-name=RASQUAL_PERM_" + self.experiment_id + "_" + str(i) + "\n")
                runsb.write("#SBATCH -o "+join(cfg.OUTPUT_DIR, "logs", self.experiment_id, "RASQUAL_PERM_" + str(i) + "_" + "%J.o")+"\n")
                runsb.write("#SBATCH -e "+join(cfg.OUTPUT_DIR, "logs", self.experiment_id, "RASQUAL_PERM_" + str(i) + "_" + "%J.e")+"\n")
                runsb.write("\n" + " ".join(command) + "\n")
            os.system("sbatch " + run_file)

class Check_All_RASQUAL_permutations(luigi.Task):

    experiment_id = luigi.Parameter()

    def run(self):

        tests = Create_RASQUAL_Tests(self.experiment_id)
        yield tests

        with open(tests.output().path, "r") as t:
            all_tests = t.readlines()

        for i,test in enumerate(all_tests):
            yield Run_RASQUAL_permutation(self.experiment_id,str(i))
