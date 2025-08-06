import argparse
import datetime
import glob
import os
import re
import subprocess


def load_arguments() -> argparse.Namespace:
    """
    Loads and parses arguments.
    :return: args - parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Python script to run the whole Dante REMASTR pipeline for several files')

    # required
    parser.add_argument('bam_files', help='Glob path to bam file(s) that needs to be processed')
    parser.add_argument('nomenclature_file', help='Path to nomenclature file')

    # additional args
    parser.add_argument('output_dir', help='Path to directory where the output will be stored. Default=<current_dir>/results', nargs='?',
                        default=None)
    parser.add_argument('--skipbams', type=int, help='Skip some BAMs. Default=0', default=0)
    parser.add_argument('--maxbams', type=int, help='Process maximally this number of BAMs. Default=all', default=None)
    parser.add_argument('--dante-repo', help='Path to Dante repository. Default=use this location', default=None)
    parser.add_argument('--reference-fasta', help='Path to reference fasta. Default=/data/genome/human/grch38_decoy/grch38_decoy.fa',
                        default='/data/genome/human/grch38_decoy/grch38_decoy.fa')
    parser.add_argument('--param-file', type=str, help='Parameter file for inference of alleles. Default=MiSeq parameters from repo',
                        default=None)
    parser.add_argument('--skip-subpath', action='store_true', help='Do not append sub-path to output as <nomenclature_name>/<sample_name>')
    parser.add_argument('--dante-args', '-a', help='Send these additional arguments to the dante_remastr call. Default=None', default='')
    parser.add_argument('--remastr-args', '-r', help='Send these additional arguments to the remastr call. Default=None', default='')
    parser.add_argument('--use-old-defaults', '-d', action='store_true', help='Use old defaults for dante_remastr call')
    parser.add_argument('--visualize', '-v', action='store_true', help='Print all the outputs. Default is to print only the result table to stdout.')
    parser.add_argument('--skip-call', '-s', action='store_true', help='Skip the calling of scripts, just do a dry run.')
    parser.add_argument('--skip-remastr', action='store_true', help='Skip the calling of remastr, use only when recomputing results.')
    parser.add_argument('--ignore-errors', '-e', action='store_true', help='Ignore errors and continue, otherwise end after first non 0 return code.')
    parser.add_argument('--additional-bams', help='More glob paths to bam file(s) to process.', nargs='*', default=[])
    parser.add_argument('--regex-bam', help='Regex for basename of bam files for finer control. Default=\'^[^.]+\\.[^.]+$\' (all with one dot)',
                        default=r'^[^.]+\.[^.]+$')
    parser.add_argument('--gzipped', action='store_true', help='Work with gzipped remastr output.')

    # preprocessing
    parser.add_argument('--flank-len', '-f', type=int, help='Flank length to include into HMM (remastr). Default=30', default=30)
    parser.add_argument('--min-mapq', '-m', type=int, help='Minimal mapping quality to accept. Use 0 for no filter. Default=30', default=30)
    parser.add_argument('--include-duplicates', action='store_true', help='Include marked duplicates. Duplication marking should be done before '
                                                                          'the pipeline and it is not handled by the pipeline.')
    parser.add_argument('--males', help='Male samples. Specify either sample names or sample serial numbers (1-based). Default=[]',
                        nargs='+', default=[])

    # converting bam file
    parser.add_argument('--bam-conversion', help='Task to convert bam file before processing. for example stripping weird headers: '
                                                 '"samtools reheader -c \'grep -v ^@RG\' {input} > {output}"', default=None)
    parser.add_argument('--bam-temp-loc', help='Temporary location for bam files. Default=<home_dir>', default='~')
    parser.add_argument('--additional-nomenclatures', help='More nomenclature files to process on the same BAM file.', nargs='*', default=[])

    args = parser.parse_args()

    if args.dante_repo is None:
        args.dante_repo = os.path.dirname(os.path.realpath(__file__))
    if args.output_dir is None:
        args.output_dir = f'{os.getcwd()}/results'
    if args.param_file is None:
        args.param_file = f'{args.dante_repo}/src/inference/training_input/miseq_params.params'
    if args.maxbams is None:
        args.maxbams = 10000000

    return args


def run_command(command: str, exit_on_error: bool, skip_call: bool) -> None:
    """
    Prints and runs the command
    :param command: str - command to run
    :param exit_on_error: bool - whether to exit on error
    :param skip_call: bool - dry run? do not call the command
    """
    start_time = datetime.datetime.now()
    print(f'#Start of a command: {start_time}')
    print(command)
    if not skip_call:
        return_code = subprocess.call(command, shell=True)
        end_time = datetime.datetime.now()
        print(f'#End   of a command: {end_time}')
        print(f'#Time elapsed      : {end_time - start_time}')
        if exit_on_error and return_code != 0:
            exit(return_code)


if __name__ == '__main__':
    # load arguments
    args = load_arguments()

    # create the output dir if needed
    os.makedirs(args.output_dir, exist_ok=True)

    # select bam files
    bam_paths = [args.bam_files] + args.additional_bams
    bam_files = []
    for bam_path in bam_paths:
        bam_files.extend(glob.glob(bam_path))
    bam_files = [bam for bam in bam_files if re.match(args.regex_bam, os.path.basename(bam))]
    bam_files = sorted(bam_files)[args.skipbams:args.skipbams + args.maxbams]

    print(f'Working on {len(bam_files)} samples...')
    print('conda activate remastr')

    # gather nomenclatures
    nomenclatures = [args.nomenclature_file] + args.additional_nomenclatures
    nomenclature_names = [os.path.basename(nomenclature).split('.')[0] for nomenclature in nomenclatures]

    # go through the bam files
    for i, bam_file in enumerate(bam_files):
        # build remastr call
        sample_name = os.path.basename(bam_file).split('.')[0]
        male_sample = str(i + 1) in args.males or sample_name in args.males

        print(f'\nSample {sample_name} ({i+1}/{len(bam_files)})')

        # recompute bam file:
        temp_bam = None
        recompute_call = None
        if args.bam_conversion is not None:
            temp_bam = f'{os.path.expanduser(args.bam_temp_loc)}/{sample_name}_temp_recomputed.bam'
            recompute_call = args.bam_conversion.format(input=bam_file, output=temp_bam)
            bam_file = temp_bam
            if not os.path.exists(temp_bam):
                run_command(recompute_call, not args.ignore_errors, args.skip_call or args.skip_remastr)

        # go through all nomenclature names
        for nomenclature_file, nomenclature_name in zip(nomenclatures, nomenclature_names):
            output_dir = args.output_dir if args.skip_subpath else f'{args.output_dir}/{nomenclature_name}/{sample_name}'
            output_file = f'{output_dir}/{sample_name}.tsv'

            # build remastr call
            remastr_params = '--flank 30 --quality 0' if args.use_old_defaults else (f'--flank {args.flank_len} --quality {args.min_mapq}' +
                                                                                     (' -d' if not args.include_duplicates else ''))
            remastr_params += ' ' + args.remastr_args
            remastr_call = (f'{args.dante_repo}/remastr/target/release/dante_cli -f {args.reference_fasta} -m {nomenclature_file} '
                            f'-b {bam_file} -o {output_file} {remastr_params}')

            # gzip results?
            gzip_call = ''
            if args.gzipped:
                gzip_call = f'gzip {output_file}'
                output_file += '.gz'

            # build dante_remastr call (report)
            additional_params = ' --processes 10 --min-flank-len 5 --min-rep-len 6 --min-rep-cnt 2 --max-rel-error 0.05' if args.use_old_defaults else ''
            additional_params += (' ' + args.dante_args + (' -v' if args.visualize else '') + (' --male' if male_sample else '') +
                                  (' --input-gzipped' if args.gzipped else ''))
            dante_call = (f'python {args.dante_repo}/dante_remastr.py --param-file {args.param_file}{additional_params} -o {output_dir} < '
                          f'{output_file} > {output_dir}/{sample_name}_res.tsv')

            # print and call commands
            if not args.skip_call:
                os.makedirs(output_dir, exist_ok=True)
            run_command(remastr_call, not args.ignore_errors, args.skip_call or args.skip_remastr)
            if args.gzipped:
                run_command(gzip_call, not args.ignore_errors, args.skip_call or args.skip_remastr)
            run_command(dante_call, not args.ignore_errors, args.skip_call)

        # remove temporary bam
        if temp_bam is not None:
            run_command(f'rm {temp_bam}*', not args.ignore_errors, args.skip_call or args.skip_remastr)

    # gather results
    for nomenclature_name in nomenclature_names:
        output_dir = args.output_dir if args.skip_subpath else f'{args.output_dir}/{nomenclature_name}'
        call_gather = f'python {args.dante_repo}/gather_result_files.py {output_dir}'
        call_motif = f'python {args.dante_repo}/create_motif_report.py {output_dir}'
        call_visual = f'python {args.dante_repo}/create_visual_report.py {output_dir}'
        run_command(call_gather, not args.ignore_errors, args.skip_call)
        run_command(call_motif, not args.ignore_errors, args.skip_call)
        run_command(call_visual, not args.ignore_errors, args.skip_call)
