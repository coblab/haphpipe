from __future__ import print_function
import os
import sys
import argparse
import shutil
import json

from past.utils import old_div
from Bio import SeqIO

from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import MissingRequiredArgument

__author__ = 'Margaret C. Steiner, Keylie M. Gibson, and Matthew L. Bendall'
__copyright__ = 'Copyright (C) 2020 Margaret C. Steiner and (C) 2019 Keylie M. Gibson and Matthew L. Bendall'

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--fq1',type=sysutils.existing_file,
                        help='Fastq file with read 1')
    group1.add_argument('--fq2', type=sysutils.existing_file,
                        help='Fastq file with read 2')
    group1.add_argument('--fqU',type=sysutils.existing_file,help='Fastq file with unpaired reads')
    group1.add_argument('--ref_fa',type=sysutils.existing_file,help="Reference FASTA file")
    group1.add_argument('--outdir',type=sysutils.existing_dir,default='.',help='Output directory')

    group2 = parser.add_argument_group('CliqueSNV Options')
    group2.add_argument('--jardir',type=str,help='Path to clique-snv.jar (existing) (Default: current directory)')
    group2.add_argument('--O22min',type=float,help="minimum threshold for O22 value")
    group2.add_argument('--O22minfreq',type=float,help="minimum threshold for O22 frequency relative to read coverage")
    group2.add_argument('--printlog',action='store_true',help="Print log data to console")
    group2.add_argument('--merging',type=str,help='Cliques merging algorithm: accurate or fast')
    group2.add_argument('--fasta_format',type=str,default='extended4',help="Fasta defline format: short or extended, add number at end to adjust precision of frequency")
    group2.add_argument('--outputstart',type=int,help="Output start position")
    group2.add_argument('--outputend', type=int, help="Output end position")

    group3 = parser.add_argument_group('HAPHPIPE Options')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                        (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Name for log file (output)')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    group3.add_argument('--ncpu', type=int, default=1, help='Number of CPU to use')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')

    parser.set_defaults(func=cliquesnv)

def cliquesnv(fq1=None,fq2=None,fqU=None,ref_fa=None,outdir='.',jardir=None,O22min=None,O22minfreq=None,printlog=None,single=False,
              merging=None,fasta_format='extended4',outputstart=None,outputend=None,keep_tmp=False, quiet=False, logfile=None, debug=False,ncpu=1):

    # check if paired vs. single
    if fq1 is None and fq2 is None and fqU is not None:
        single=True

    # check dependencies and required arguments
    if fq1 is None and fq2 is None and fqU is None:
        raise MissingRequiredArgument("No fastq files given.")
    if single == False and (fq1 is None or fq2 is None):
        raise MissingRequiredArgument("Either fq1 or fq2 missing.")
    if ref_fa is None:
        raise MissingRequiredArgument("Reference FASTA missing.")

    sysutils.check_dependency('samtools')
    sysutils.check_dependency('bwa')

    if jardir is None:
        # jardir was not provided, assume cliquesnv installed from conda
        USE_JAR = False        
        sysutils.check_dependency('cliquesnv')
    else:
        USE_JAR = True
        if(os.path.isfile(os.path.join(jardir,"clique-snv.jar"))):
            print("CliqueSNV JAR file found.")
        else:
            raise MissingRequiredArgument("No JAR file found.")

    # Temporary directory
    tempdir = sysutils.create_tempdir('clique_snv', None, quiet, logfile)

    # Load reference fasta
    refs = {s.id: s for s in SeqIO.parse(ref_fa, 'fasta')}

    # Identify reconstruction regions
    regions = []
    for rname, s in refs.items():
        regions.append(('cs%02d' % (len(regions) + 1), rname, 1, len(s)))

    sysutils.log_message('[--- Haplotype Reconstruction Regions ---]\n', quiet, logfile)
    for iv in regions:
        sysutils.log_message('%s -- %s:%d-%d\n' % iv, quiet, logfile)

    if single == False: #paired end
        # remove .1 and .2 from read names
        fq1_c = os.path.join(tempdir,"fq1_corrected.fastq")
        fq2_c = os.path.join(tempdir, "fq2_corrected.fastq")
        cmd01 = ["cat %s | sed 's/\.1 / /' > %s" % (fq1,fq1_c)]
        cmd02 = ["cat %s | sed 's/\.2 / /' > %s" % (fq2,fq2_c)]
        sysutils.command_runner([cmd01,cmd02],'clique_snv:setup',quiet,logfile,debug)

        # Create alignment for each REFERENCE in the reconstruction regions
        alnmap = {}
        for cs, rname, spos, epos in regions:
            if rname not in alnmap:
                # Create alignment
                tmp_ref_fa = os.path.join(tempdir, 'ref.%d.fa' % len(alnmap))
                tmp_sam = os.path.join(tempdir, 'aligned.%d.sam' % len(alnmap))
                SeqIO.write(refs[rname], tmp_ref_fa, 'fasta')
                cmd1 = ['bwa', 'index', tmp_ref_fa, ]
                cmd2 = ['bwa', 'mem', tmp_ref_fa, fq1_c, fq2_c, '|', 'samtools', 'view', '-h', '-F', '12', '>', tmp_sam, ]
                cmd3 = ['rm', '-f', '%s.*' % tmp_ref_fa]
                sysutils.command_runner(
                    [cmd1, cmd2, cmd3], 'clique_snv:setup', quiet, logfile, debug
                )
                alnmap[rname] = (tmp_ref_fa, tmp_sam)

    else: #single read

        # Create alignment for each REFERENCE in the reconstruction regions
        alnmap = {}
        for cs, rname, spos, epos in regions:
            if rname not in alnmap:
                # Create alignment
                tmp_ref_fa = os.path.join(tempdir, 'ref.%d.fa' % len(alnmap))
                tmp_sam = os.path.join(tempdir, 'aligned.%d.sam' % len(alnmap))
                SeqIO.write(refs[rname], tmp_ref_fa, 'fasta')
                cmd1 = ['bwa', 'index', tmp_ref_fa, ]
                cmd2 = ['bwa', 'mem', tmp_ref_fa, fqU, '|', 'samtools', 'view', '-h', '-F', '12', '>',
                        tmp_sam, ]
                cmd3 = ['rm', '-f', '%s.*' % tmp_ref_fa]
                sysutils.command_runner(
                    [cmd1, cmd2, cmd3], 'clique_snv:setup', quiet, logfile, debug
                )
                alnmap[rname] = (tmp_ref_fa, tmp_sam)


    # Run CliqueSNV for each region
    cmd4 = ['mkdir -p %s' % os.path.join(outdir, 'clique_snv')]
    sysutils.command_runner([cmd4, ], stage='cliquesnv', quiet=quiet, logfile=logfile, debug=debug)
    i=0 #index for filenames
    for cs, rname, spos, epos in regions:
        msg = "Reconstruction region %s:" % cs
        msg += " %s:%d-%d\n" % (rname, spos, epos)
        sysutils.log_message(msg, quiet, logfile)

        # rename the cliquesnv number (cs##) to include region (now: cs##_reg)
        cs = '%s_%s' % (cs, rname.split('|')[-2])

        samfile = os.path.join(tempdir,'aligned.%d.sam' % i)
        method = 'snv-illumina'
        if USE_JAR:
            cmd5 = ['java', '-jar', os.path.join(jardir,'clique-snv.jar'), ]
        else:
            cmd5 = ['cliquesnv', ]
        
        cmd5 += [
            '-m %s' % method,
            '-in %s' % samfile,
            '-threads %d' % ncpu,
            '-outDir %s' % tempdir,
            '-fdf %s' % fasta_format,
        ]
        if O22min is not None:
            cmd5 += ['-t %f' % O22min]
        if O22minfreq is not None:
            cmd5 += ['-tf %f' % O22minfreq]
        if printlog is not None:
            cmd5 += ['-log']
        if merging is not None:
            cmd5 += ['-cm %s' % merging]
        if outputstart is not None:
            cmd5 += ['-os %d' % outputstart]
        if outputend is not None:
            cmd5 += ['-oe %d' % outputend]
        sysutils.command_runner([cmd5, ], stage='clique_snv', quiet=quiet, logfile=logfile, debug=debug)

        # copy output file and delete tempdir
        os.makedirs(os.path.join(outdir, 'clique_snv/%s' % cs), exist_ok=True)
        
        out_txt = os.path.join(outdir,'clique_snv/%s/%s.txt' % (cs,cs))
        out_fas = os.path.join(outdir,'clique_snv/%s/%s.fasta' % (cs,cs))
        out_json = os.path.join(outdir,'clique_snv/%s/%s.json' % (cs,cs))
        out_sum = os.path.join(outdir,'clique_snv/%s/%s_summary.txt' % (cs,cs))
        
        if os.path.exists(os.path.join(tempdir, 'aligned.%d.txt' % i)):
            shutil.copy(os.path.join(tempdir, 'aligned.%d.txt' % i), out_txt)
        
        if os.path.exists(os.path.join(tempdir, 'aligned.%d.fasta' % i)):
            shutil.copy(os.path.join(tempdir, 'aligned.%d.fasta' % i), out_fas)
        
        if os.path.exists(os.path.join(tempdir, 'aligned.%d.json' % i)):
            shutil.copy(os.path.join(tempdir, 'aligned.%d.json' % i), out_json)

        # parse output file, may be txt or json
        if os.path.exists(out_txt):
            with open(out_sum, 'w') as sumfile, open(out_txt, 'r') as infile:
                l = infile.readlines()
                freqs = []
                haps = []
                tempnum=''
                for line in l:
                    if "SNV got" in line:
                        tempnum = line.split(' ')[2]
                    if "frequency" in line:
                        freqs += [float(line.split(' ')[2][:-2])]
                    if "haplotype=" in line:
                        haps += [line.split('=')[1][1:-2]]
                sumfile.write('CliqueSNV_num_hap\t%s\n' % tempnum)
    
                freq_sqrd = [x ** 2 for x in freqs]
                freq_sqrd_sum = sum(freq_sqrd)
                hap_div = ((old_div(7000, (7000 - 1))) * (1 - freq_sqrd_sum))
                sumfile.write('CliqueSNV_hap_diversity\t%s\n' % hap_div)
                sumfile.write('CliqueSNV_seq_len\t%s\n' % len(haps[0]))
        elif os.path.exists(out_json):
            with open(out_sum, 'w') as sumfile, open(out_json, 'r') as infile:
                dat = json.load(infile)
                sumfile.write('CliqueSNV_num_hap\t%s\n' % dat['foundHaplotypes'])
                freqs = [h['frequency'] for h in dat['haplotypes']]
                freq_sqrd = [x ** 2 for x in freqs]
                freq_sqrd_sum = sum(freq_sqrd)
                hap_div = ((old_div(7000, (7000 - 1))) * (1 - freq_sqrd_sum))           
                sumfile.write('CliqueSNV_hap_diversity\t%f\n' % hap_div)
                sumfile.write('CliqueSNV_seq_len\t%s\n' % len(dat['haplotypes'][0]['haplotype']))        
        
        with open(os.path.join(outdir, 'clique_snv/%s/%s.fasta' % (cs, cs)), 'r') as fastafile:
            fastadata=fastafile.read().replace('aligned.%d' % i,rname)
        
        with open(os.path.join(outdir, 'clique_snv/%s/%s.fasta' % (cs, cs)), 'w') as newfastafile:
            newfastafile.write(fastadata)

        i += 1

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'clique_snv', quiet, logfile)

    return

def console():
    parser = argparse.ArgumentParser(
        description='Haplotype reconstruction with CliqueSNV.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    try:
        args.func(**sysutils.args_params(args))
    except MissingRequiredArgument as e:
        parser.print_usage()
        print('error: %s' % e, file=sys.stderr)


if __name__ == '__main__':
    console()
