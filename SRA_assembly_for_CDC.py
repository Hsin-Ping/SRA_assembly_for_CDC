import time
import os
import re
import sys
import shutil
import argparse
import subprocess
import xml.etree.cElementTree as ET
from tempfile import TemporaryDirectory
import pandas as pd
import time
import glob

def run_cmd(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, check=True)
    return p

def run_cmd2(cmd):
    #p = subprocess.run(cmd, stderr=subprocess.PIPE, shell=True, check=True)
    p = subprocess.run(cmd, shell=True, check=True)
    return p

###for_assembly
def bases_percentage(filepath, qscore=0):
    p = run_cmd(f"seqtk fqchk -q {qscore} {filepath} | grep ALL | awk '{{print $NF}}'")
    return float(p.stdout)

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
ADAPTERS = os.path.join(CURRENT_DIR, 'trimmomatic.fa')
MIN_BQ = 3

def crop_position(filepath, window_size=3, gap=10):
    p = run_cmd(f"seqtk fqchk {filepath}")
    fq_check = p.stdout.decode().strip().split('\n')[3:]
    fq_check = (line.split()[2:6] for line in fq_check)
    content_gaps = []
    for line in fq_check:
        a, c, g, t = [float(value) for value in line]
        content_gaps.append(max(abs(a - t), abs(c - g)))
    # check from forward
    for start in range(len(content_gaps)):
        end = start + window_size
        window = content_gaps[start: end]
        if max(window) < gap:
            headcrop = start
            break
        else:
            headcrop = 0
    # check from revers
    for start in range(len(content_gaps), 0, -1):
        end = start - window_size
        window = content_gaps[end:start]
        if max(window) < 10:
            crop = start
            break
        else:
            crop = len(content_gaps)
    return crop, headcrop

def trimming(forward_reads, reverse_reads, outdir, threads):
    crop, headcrop = crop_position(forward_reads)
    opt = f"CROP:{crop} HEADCROP:{headcrop} ILLUMINACLIP:{ADAPTERS}:2:30:10 LEADING:{MIN_BQ} TRAILING:{MIN_BQ} " \
          f"SLIDINGWINDOW:4:20 MINLEN:36 TOPHRED33"
    paired_1 = os.path.join(outdir, 'R1.fq')
    paired_2 = os.path.join(outdir, 'R2.fq')
    cmd = f"java -jar trimmomatic-0.39.jar PE -threads {threads} {forward_reads} {reverse_reads} {paired_1} /dev/null" \
          f" {paired_2} /dev/null {opt}"
    run_cmd(cmd)
    return paired_1, paired_2

def dump_fastq_from_sra(srafile, outdir):
    run_cmd(f'fastq-dump --split-files --outdir {outdir} {srafile}')


class SequenceReadArchive:
    def __init__(self,filepath):
       self._set_filepath(filepath)
       self._get_stat()

    def _set_filepath(self,filepath):
        if os.access(filepath, os.F_OK) is False:
           raise FileNotFoundError("File not found.")
        with open(filepath,'rb') as handle:
           if handle.read(8).decode() == 'NCBI.sra' is False:
              raise Exception(f"File format is not 'NCBI.sra'.")
        self._filepath = filepath

    def _get_stat(self):
        p = run_cmd(f'sra-stat -x -s -b 1 -e 2 {self._filepath}')
        self._stat_tree = ET.fromstring(p.stdout.decode())

    @property
    def filepath(self):
        return self._filepath

    @property
    def layout(self):
        return self._stat_tree.find('Statistics').attrib['nreads']
###

def prefetch_sra(sralist,outdir):
    ss = " ".join(sralist)
    print("now download",ss,"runs.")
    cmd = "prefetch "+ss+" --output-directory "+outdir
    run_cmd2(cmd)

def run_dump(need_run,sra_dir,assem_dir,threads,gsize):
    k = list(range(0,len(need_run),3))
    for i in k:
        run_id = need_run[i:i+3]
        print("run_id",run_id)
        prefetch_sra(run_id,sra_dir)
        time.sleep(1)
        for x in run_id:
            start = time.time()
            path = sra_dir+"/"+x+"/*"
            re_path = "".join(glob.glob(path))
            sra_file = os.path.abspath(re_path)
            #print(x)
            #print("sra_abspath:",sra_file)
            try:
                sra = SequenceReadArchive(sra_file)
                #print("layout:",sra.layout)
            except Exception as e:
                print("error")
                #sys.exit(e)
            if sra.layout != '2':
                sys.exit(f'File layout is not pair-end')
            print('Now assembly',sra_file,'...')
            outdir = assem_dir+"/"+"".join(x)
            #print(outdir)
            fastq_dir = os.path.join(outdir,'fastq')
            os.makedirs(fastq_dir, exist_ok = True)
            print('Dump fastq.')
            dump_fastq_from_sra(sra_file, fastq_dir)
            forward_reads, reverse_reads = [os.path.join(fastq_dir,i) for i in os.listdir(fastq_dir)]
            print('Trim sequences.')
            r1, r2 = trimming(forward_reads, reverse_reads, fastq_dir, threads)
            if bases_percentage(r1, 30) < 90 and bases_percentage(r2, 30) < 90:
                shutil.rmtree(outdir)
                sys.exit('Reads quality is too low.')
            print("Run assembly pipline 'shovill'")
            assemble_dir = os.path.join(outdir,'assembly_result')
            cmd = f"shovill --R1 {r1} --R2 {r2} --outdir {assemble_dir} --depth 100 --tmpdir . --cpus {threads} --ram 3 --force"
            if gsize:
                cmd += f" --gsize {gsize}"
            #print(cmd)
            run_cmd(cmd)
            print('Done,total cost',time.time()-start,'secs')
        shutil.rmtree(sra_dir)
        os.mkdir(sra_dir)

def main():
    parser = argparse.ArgumentParser("Download_avaliable_runs_from_NCBI_sra_database.")
    parser.add_argument("--avaliable",required = True, help = "")
    parser.add_argument("--downloaded",required = True, help = "")
    parser.add_argument("--sra_dir",required = True, help ="temp folder to downlaod runs.")
    parser.add_argument("--assembly_dir", required = True, help = "folder to save assembly")
    parser.add_argument("--tmpdir", default="/tmp", help="Directory of temp folder default: '/tmp'")
    parser.add_argument("--gsize", default='', help="Estimated genome size(MB) eg. 3.2M. default: ''")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    args = parser.parse_args()
    
    ava = args.avaliable
    don = args.downloaded
    sra_dir = args.sra_dir
    assem_dir = args.assembly_dir
    tmpdir = args.tmpdir
    threads = args.threads
    gsize = args.gsize

    ava_df = pd.read_csv(ava)
    ava_run = list(ava_df['Run'])
    don_df = pd.read_table(don)
    don_run = list(don_df['Run'])
    need_run = list(filter(lambda x : x not in don_run,ava_run))
    print("Toal",len(need_run),"sra runs need to downlaod.")
    os.makedirs(sra_dir, exist_ok = True)
    os.makedirs(assem_dir,exist_ok = True)
    #sra_dir = os.path.abspath(sra_dir)
    run_dump(need_run,sra_dir,assem_dir,threads,gsize)
    
if __name__ == '__main__':
   main()
