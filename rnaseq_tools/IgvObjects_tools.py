import os

def run_IGV_script(igv_script, igv_jar, memMB):
    '''
    Run an IGV batch script
    '''
    import datetime
    # get the X11 Xvfb port number
    x_serv_port = get_open_X_server()
    # build the system command to run IGV
    # igv_command = "(Xvfb :{} &) && DISPLAY=:{} java -Xmx{}m -jar {} -b {} && killall Xvfb".format(x_serv_port, x_serv_port, memMB, igv_jar, igv_script)
    igv_command = "xvfb-run --auto-servernum --server-num=1 java -Xmx{}m -jar {} -b {}".format(memMB, igv_jar, igv_script)
    os.system(igv_command) # not with subcommand -- no need to make script wait until finished


def takeSnapshot(bam_file, region_file, IGV_batchfile_script, output_dir, genome,
                 image_height = 500, nf4_mode = True, fig_format = 'png'):
    """
    Downloaded from: https://github.com/stevekm/IGV-snapshot-automator
    Modified by Yiming Kang
    Revised by chase mateusiak

    Per the README on the github above:

    This script will load IGV in a virtual X window, load all supplied input files
    as tracks, and take snapshots at the coorindates listed in the BED formatted
    region file.

    example IGV batch script:

    new
    snapshotDirectory IGV_Snapshots
    load test_alignments.bam
    genome hg19
    maxPanelHeight 500
    goto chr1:713167-714758
    collapse
    snapshot chr1_713167_714758_h500.[png/svg]
    goto chr1:713500-714900
    collapse
    snapshot chr1_713500_714900_h500.[png/svg]
    exit
    """

    # write the IGV batch script
    write_IGV_script(input_files=input_files, region_file=region_file, IGV_batchscript_file=batchscript_file,
                     IGV_snapshot_dir=output_dir, genome_version=genome, image_height=image_height,
                     nf4_mode=nf4_mode, fig_format=fig_format)


    run_IGV_script(igv_script=batchscript_file, igv_jar=igv_jar_bin, memMB=igv_mem)


parser.add_argument("-ht", default='500', type=str, dest='image_height', metavar='image height',
                    help="Height for the IGV tracks")
parser.add_argument("-o", default='IGV_Snapshots', type=str, dest='outdir', metavar='output directory',
                    help="Output directory for snapshots")
parser.add_argument("-bin", default="bin/IGV_2.3.81/igv.jar", type=str, dest='igv_jar_bin', metavar='IGV bin path',
                    help="Path to the IGV jar binary to run")
parser.add_argument("-mem", default="4000", type=str, dest='igv_mem', metavar='IGV memory (MB)',
                    help="Amount of memory to allocate to IGV, in Megabytes (MB)")
parser.add_argument("-nosnap", default=False, action='store_true', dest='no_snap',
                    help="Don't make snapshots, only write batchscript and exit")
parser.add_argument("-suffix", default=None, dest='suffix',
                    help="Filename suffix to place before '.png' in the snapshots")
parser.add_argument("-nf4", default=False, action='store_true', dest='nf4_mode',
                    help="'Name field 4' mode; uses the value in the fourth field of the regions file as the filename for each region snapshot")
parser.add_argument("-onlysnap", default=False, dest='onlysnap',
                    help="Path to batchscript file to run in IGV. Performs no error checking or other input evaluation, only runs IGV on the batchscript and exits.")
parser.add_argument("-fig_format", default='png', dest='fig_format',
                    help="Image format of snapshot. Choose from png (default) or svg.")

args = parser.parse_args()

input_files = args.input_files
region_file = args.region_file
genome = args.genome
image_height = args.image_height
outdir = args.outdir
igv_jar_bin = args.igv_jar_bin
igv_mem = args.igv_mem
no_snap = args.no_snap
suffix = args.suffix
nf4_mode = args.nf4_mode
onlysnap = args.onlysnap
fig_format = args.fig_format

main(input_files=input_files, region_file=region_file, genome=genome, image_height=image_height,
     outdir=outdir, igv_jar_bin=igv_jar_bin, igv_mem=igv_mem, no_snap=no_snap, suffix=suffix,
     nf4_mode=nf4_mode, onlysnap=onlysnap, fig_format=fig_format)

