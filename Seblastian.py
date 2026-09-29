#! /usr/bin/env python2.7
__author__  = "Marco Mariotti"
__email__   = "marco.mariotti@ub.edu"
__version__ = "1.2a"
from string import *
import sys, os
import traceback
#sys.path.insert(0, "/users-d3/mmariotti/libraries/selenoprofiles_2/")
#sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import *
#from profiles_classes import exonerated, partial_to_profile_alignment

from selenoprofiles_2 import blasthit, load_chromosome_lengths, exoneratehit, superexoneratehit, set_selenoprofiles_var, get_selenoprofiles_var, parse_secisearch, merge_p2g_hits_by_colinearity, secisearchhit, cyclic_exonerate

help_msg="""This program is a method for identification of selenoproteins in genomes. It initially searches the genome for SECIS elements, then it looks upstream for sequences likely to be coding for selenoproteins: it performs a blastx search against a protein database and looks for alignments which look like selenoproteins upstream of the SECIS. It finally runs exonerate to define the whole gene structures of the candidates.

Usage:   Seblastian.py   target_file.fasta   output_prefix   [ protein_database.fasta ]

### Routines:     these are run by default if the corresponding output file is not present already. Use these options to force their execution.
-search     || -S   run methods for searching secis on the target file
-extract    || -X   extract the sequence upstream of each secis
-find_uga   || -U   annotate the position of each UGA in the upstream sequences
-blast      || -B   run blastx, the query is all upstream sequences (only positive strand is translated), the target is the protein database. The blast output is also parsed
-fetch      || -T   fetch all protein sequences which appear in the filtered blast output and save the gi codes for those containing U (known selenoproteins)
-parse      || -P   parse the blast output looking for TGA to Sec/Cys/XXX alignments or for alignments with any annotated selenoprotein or for other proteins ending just upstream
-exonerate  || -E   run exonerate to find the complete gene structures of the candidates
-filter     || -F   choose either the exonerate or reverseblastx prediction and filter them
-output     || -O   generate output files 

### Options per category (+ indicates an argument is required):
  ## search 
-no_complement  || -c     to turn off the search in the complementary strand (for RNA data)
-method  || -m   +[,+]    method to predict secis. Multiple comma separated argument are accepted. possibilities: infernal (default), secisearch1, covels
-no_secis_filter || -nf   filter secis using the custom built-in filter.
-secis_energy  +          energy threshold for secis filtering (used only if secis filter is active)
-no_refine_ss             do not refine secondary structure of SECIS_covels
-type                     try to predict stem3, to discriminate between type I and II SECIS. This information will then be displayed in the summary of the SECIS
-load_secis               force loading of secis from the output files and not from dumped, processed secis file. Useful if you want to change filtering of SECIS but not rerun SECISearch

  # for infernal method
-infernal_score  +        use a score threshold for infernal
-infernal_evalue +        use a evalue threhold for infernal
-infernal_cm     +        path to the infernal secis model searched
-infernal_stk    +        path to the stk alignment used to build the secis model searched
-n_secis +                keep only the best n secis found in the search output (after filtering, if active)
-infernal_no_mpi || -i    do not try to use mpi to run cmsearch with multiple processors
  # for secisearch method
-secisearch_pattern || -p +   use this as pattern for patscan phase of SECISearch1. possible values: default, strict, loose, loose+
  # for covels method
-covels_cm        +       path to the covels secis model searched
-covels_score     +       use a score threshold during covels search
-covels_extension +       the extension done on each side of the cove predicted range before aligning to the infernal model

  ## extract    
-length_gap  +            length of the sequence just upstream of the secis which is NOT fetched (to forbid the coding sequence being adjacent to the SECIS) [def:0]
-length_upstream  +       length of the upstream sequence which is fetched and then blastxed [def:5000]

  ## blast      
-net  || -n               use netblast instead of local blastall (blastcl3 must be installed) ### not working yet!!
-blast_evalue +           evalue threshold used with blastall [def: 1e-3]
-blast_opt "+"            additional options to be used with blastall
-blast_filter  "+"        further filter blast hits using a python procedure. you can refer to the blasthit as x (example: -blast_filter "x.chromosome=='chr1'")

  ## parse
-max_unaligned_end  +     maximum aminoacids unaligned at the c terminal to consider a blast hit for labelling it as a CDS ending just upstream

  ## exonerate
-cyclic_exonerate_ext +   extension used for the cyclic exonerate routine (see selenoprofiles)
-keep_blast_labels  +,+   labels for which the exonerate hit will be ignored, and the blast information will be passed to the filter instead

  ## filter
-max_secis_distance +     set the maximum CDS-secis distance for results to pass the filter
-output_labels     +,+    only results with any of these labels will be output. the order of these labels determine the order in output
-custom_filter     "+"    filter results using a defined python procedure instead of the default filtering procedure. You can refer to the result object as x (example: -blast_filter "x.chromosome=='chr1'")                

  ## output 
if the following options are active, the corresponding output format is produced. the default value is shown inside brackets. Use 0 as argument to turn off an option
-out_pretty      [1]      a text file with a summary including the protein prediction and the secis structure for all filtered results
-out_gff         [1]      a text GFF file with the coordinates of all predicted proteins and their secis for all filtered results
-out_fasta       [0]      a text FASTA file with the predicted protein sequence for all filtered results
-out_cds         [0]      a text FASTA file with the predicted coding sequence for all filtered results
-out_secis_fasta [0]      a text FASTA file with the predicted nucleotide sequence for the SECIS elements of all filtered results 
-out_secis_image || -I   [0]      a PNG image for each SECIS elements of the filtered results is produced in a dedicated folder
-fd                       SECIS image faces down instead of up
-out_all_secis   [0]      write a output file with all predicted SECIS elements, regardless of whether a selenoprotein was be found upstream of them
-out_filtered_secis [0]   write a output file with all filtered out SECIS elements, including a line with the reason why it was filtered out
-dpi_resolution  +        define the resolution for SECIS images in dpi (dot per inches). Default: 150
-three_prime_length +     max length of three prime considered when scanning for next stop codon. this is just printed in "pretty" output

  ## system and miscellaneous
-genetic_code +           NCBI identifier for genetic code (e.g. 1=standard, 2=mitochondrial ...)
-exclude      +           file with SECIS predictions to exclude. Accepted formats: .gff, .infernal, .covels, .secisearch
-processors   +           number of processors, used with infernal and blastall 
-temp  +                  temporary folder. A folder with random name is created here, used and deleted at the end of the computation.
-save_chromosomes         if active, each chromosome file which is fetched is stored in the temporary folder, in a folder named after the target file (active by default, run with "0" as argument to turn off)
-print_opt                print currently active options
-time                     print time statistics on the various phases
-bin_folder  +            looks for executables in this folder (adding it to $PATH before executing any command)
-print_commands           print all bash commands to screen before running
-v                        verbose mode
-only_search || -SS       use program only for secis prediction: it stops after searching for them with active methods (activates option -out_all_secis)
-h OR --help              print this help and exit

"""
command_line_synonyms={'S':'search', 'X':'extract', 'U':'find_uga', 'B':'blast', 'T':'fetch', 'P':'parse', 'E':'exonerate', 'O':'output', 'p':'secisearch_pattern', 'm':'method', 'n':'net', 'F':'filter', 'c':'no_complement', 'i':'infernal_no_mpi', 'e':'check_middle', 'nf':'no_secis_filter', 'SS':'only_search', 'I':'out_secis_image'}

### default options
def_opt= {'temp':'/opt/~/Seblastian/temp/', 
't':0, 'd':'/opt/~/Seblastian/uniref50.only_selenoproteins.fa', 
'o':0, 
'v':0, 'Q':0,
'search':0,   'find_uga':0, 'extract':0, 'blast':0, 'fetch':0, 'parse':0, 'exonerate':0,  'filter':0, 'output':0, 'load_secis':0, 
'method':'infernal',
'infernal_score':10, 'infernal_evalue':None, 'infernal_cm':'/opt/~/Seblastian/SECIS_infernal.cm', 'infernal_stk':'/opt/~/Seblastian/SECIS_infernal.stk',
'infernal_no_mpi':1,
'n_secis':0, 
'secisearch_pattern':'default',
'covels_cm':'/opt/~/Seblastian/SECIS_covels.cm',
'covels_score':5,
'no_secis_filter':0, 'secis_energy':-4, 'no_refine_ss':0, 'type':0,
'processors':6,
'no_complement':0, 
'save_chromosomes':1,
'length_upstream':5000, 'length_gap':50,
'blast_evalue':'1e-3', 'blast_filter':None, 'blast_opt':'', 'net':0,
'max_unaligned_end':15, 'check_middle':0,
'cyclic_exonerate_ext':30000,
'output_labels':'knownSP,newSP',
'keep_blast_labels':'unaligned,empty,lost_sec_candidate',
'custom_filter':None, 'max_secis_distance': 3000,
'out_gff':1, 'out_fasta':0, 'out_pretty':1, 'out_secis_fasta':0, 'out_all_secis':0, 'out_cds':0, 'out_secis_image':0, 'out_filtered_secis':0,
'debug':0, 'print_commands':0, 'time':1,
'bin_folder':'/opt/~/Seblastian/bin',
'only_search':0,
'color_mode':2,
'three_prime_length':1000,
'dpi_resolution':150, 
'fix':0, 'covels_extension':16,
'fd':0, 'genetic_code':0,
}

#########################################################
###### start main program function
category_per_pair=symmetrical_hash(); 
category_per_pair['U']={'A':'a'}
category_per_pair['G']={'A':'k', 'C':'c', 'U':'w'}  
category_per_pair_hierarchy='cawko'
all_methods=['infernal', 'secisearch1', 'covels', 'gff'] #, 'erpin']
executable_per_methods={'infernal':'cmsearch', 'secisearch1':'SECIScript.pl', 'covels':'covels-SE', 'gff':'echo'} #, 'erpin':'erpin'}
all_secisearch_patterns=['default', 'strict', 'loose', 'loose+', 'new']
class nohitsException(Exception):
  """ exception class to stop run if no hits are left in workflow """
def main(opt_hash={}):
#########################################################
############ loading options
  ## reading options from def_opt and from in command line, or from argument if main is run within another script
  global opt; 
  if opt_hash:  ## argument was provided
    opt=options();
    for k in def_opt: 
      if opt_hash.has_key(k): opt[k]=opt_hash[k]
      else: opt[k]=def_opt[k]
  else:  opt=command_line(def_opt, help_msg, 'tod', synonyms=command_line_synonyms, strict=1 ) # standard: running from command line
  set_MMlib_var('opt',opt)  
  #target
  global target_file;   target_file=opt['t']; check_file_presence(target_file, 'target_file');  verbose('target_file : '+target_file+'\n')
  #output
  output_prefix=opt['o']
  if not output_prefix:        
    if '.' in base_filename(target_file): output_prefix=join(target_file.split('.')[:-1], '.')
    else: output_prefix=target_file
  if output_prefix[-1]!='.':   output_prefix+='.'  
  verbose('output_prefix : '+output_prefix+'\n')  
  #setting log_file handler
  log_filename=output_prefix+'log'
  global log_file; log_file=open(log_filename, 'w'); set_MMlib_var('log_file', log_file)
  #setting print_commands 
  if opt['print_commands']: set_MMlib_var('print_commands', 1)
  #setting bin_folder if specified
  if opt['bin_folder']:     
    global bin_folder;  bin_folder=Folder(opt['bin_folder']);     set_MMlib_var('bin_folder', bin_folder)

  #input
  verbose('# INPUT FILES\n')
  global prot_database; prot_database=opt['d']; 
  #if not opt['net'] and prot_database=='nr' and bash('ls '+prot_database)[0] :    prot_database = '/users/rg/mmariotti/Databases/nr.fa' #### DEFAULT NR PATH! 
  if not opt['only_search'] and not opt['net']:  
    check_file_presence(prot_database, 'prot_database')
    # formatting db if necessary
    if bash('ls '+prot_database+'.p*')[0]: 
      service("  ...formatting database: "+prot_database)
      b=bash('formatdb -i '+prot_database+' -p T -o T')
      if b[0]:      raise Exception, "ERROR protein database cannot be formatted with formatdb: "+prot_database 
  elif opt['net']: #opt -net is active; we need to import a program to fecth sequences instead of using fastacmd to extract it from the protein database
    import fetch_ncbi_sequences

  verbose('prot_database : '+prot_database+'\n')
  # temp folders
  global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder')
  global split_folder;    
  if opt['save_chromosomes']:      split_folder=Folder(opt['temp']);               
  else:                            split_folder=temp_folder
  test_writeable_folder(split_folder); 
  set_MMlib_var('temp_folder', temp_folder); set_MMlib_var('split_folder', split_folder);
  set_selenoprofiles_var('temp_folder', temp_folder); set_selenoprofiles_var('split_folder', split_folder);    

  ## setting and printing (verbose mode) other options
  verbose('# OPTIONS'+'\n')  
  verbose('temp_folder : '+temp_folder+'\n'); verbose('split_folder : '+split_folder+'\n')
  #other options
  processors=opt['processors'];                     verbose('processors : '+str(processors)+'\n')
  blast_opt=str(opt['blast_opt']);                  verbose('blast_opt : '+str(blast_opt)+'\n')
  length_upstream=opt['length_upstream'];           verbose('length_upstream : '+str(length_upstream)+'\n')
  global length_gap; length_gap=opt['length_gap'];                     verbose('length_gap : '+str(length_gap)+'\n')
  verbose('blast_filter : '+str(opt['blast_filter'])+'\n')
  blast_evalue=opt['blast_evalue'];                 verbose('blast_evalue : '+str(blast_evalue)+'\n')
  global max_unaligned_end; max_unaligned_end=opt['max_unaligned_end'];       verbose('max_unaligned_end : '+str(max_unaligned_end)+'\n')
  cyclic_exonerate_ext=opt['cyclic_exonerate_ext']; verbose('cyclic_exonerate_ext : '+str(cyclic_exonerate_ext)+'\n')
  output_labels=opt['output_labels'].split(',');    verbose('output_labels: '+str(output_labels)+'\n')
  keep_blast_labels=opt['keep_blast_labels'].split(',');   verbose('keep_blast_labels: '+str(keep_blast_labels)+'\n')
  #output filenames
  verbose('# OUTPUT FILES'+'\n')  # # # # # # # #

  if opt['genetic_code']:       
    blast_opt+=' -Q '+str(opt['genetic_code']) 
    exonerate_version=float(bash('exonerate --version')[1].split('\n')[0].split()[-1][:3])
    if exonerate_version<2.4: 
      raise Exception, "ERROR exonerate version detected: {}\nAlternative genetic codes are bugged in exonerate versions <2.4!\nPlease download and install exonerate version 2.4.0 if you want to use -genetic_code".format(exonerate_version)

  #method
  if opt['method']=='all': opt['method']= join(  [m for m in all_methods if m!='gff'], ',')
  active_methods=opt['method'].split(',')
  for method in active_methods: 
    if not method in all_methods:  raise Exception, "ERROR method not recognized: "+method
    b=bash(executable_per_methods[method])
    if 'not found' in b[1]: raise Exception, "ERROR method not available: "+method+' ; this executable could not be found: '+executable_per_methods[method]

  #if 'not found' in bash('RNAfold -h')[1]: raise Exception, "ERROR  program RNAfold not found. You can still run with -no_refine_ss "
  if 'not found' in bash('RNAeval -h')[1]: raise Exception, "ERROR  program RNAeval not found. "

  global infernal_cm;   infernal_cm=opt['infernal_cm'];        verbose('infernal_cm : '+infernal_cm+'\n');       check_file_presence(infernal_cm, 'infernal_cm')     #model file
  global infernal_stk;  infernal_stk=opt['infernal_stk'];                 verbose('infernal_stk : '+infernal_stk+'\n');     check_file_presence(infernal_stk, 'secis_model_stk')
  global infernal_model_alignment;   infernal_model_alignment=alignment(); infernal_model_alignment.load_stockholm(infernal_stk)
  global covels_cm;    covels_cm=opt['covels_cm']; verbose('covels_cm: '+covels_cm+'\n');   check_file_presence(covels_cm, 'covels_cm')

  if 'infernal' in active_methods:
    infernal_output_file=output_prefix+'secis.infernal';                    verbose('infernal_output_file : '+infernal_output_file+'\n')
  if 'secisearch1' in active_methods:
    secisearch_pattern=opt['secisearch_pattern']; verbose('secis_pattern : '+secisearch_pattern+'\n'); 
    if not secisearch_pattern in all_secisearch_patterns: raise Exception, "ERROR secisearch pattern not recognized: "+secisearch_pattern
    secisearch_output_file=output_prefix+'secis.secisearch1';                  verbose('secisearch_output_file : '+secisearch_output_file+'\n')
  if 'covels' in active_methods:
    covels_output_file=output_prefix+'secis.covels';         verbose('covels_output_file : '+covels_output_file+ '\n')
  if 'gff' in active_methods:
    secis_gff_infile=output_prefix+'secis.gff'
  dumped_secis_filename=output_prefix+'pickled_dmp.all_secis';                verbose('dumped_secis_filename'+dumped_secis_filename+'\n')
  upstream_seqs_file=output_prefix+'upstream_seqs';                        verbose('upstream_seqs_file : '+upstream_seqs_file+'\n')
  uga_positions_file=output_prefix+'uga_pos';                              verbose('uga_positions_file : '+uga_positions_file+'\n')
  blastx_output_file=output_prefix+'blastx';                               verbose('blastx_output_file : '+blastx_output_file+'\n')
  parsed_blastx_output_file=output_prefix+'blastx.parsed';                 verbose('parsed_blastx_output_file : '+parsed_blastx_output_file+'\n')
  fetched_sequences_gis_file=output_prefix+'known_selenoprotein_gis';      verbose('fetched_sequences_gis_file : '+fetched_sequences_gis_file+'\n')
  candidates_parsed_tblastx_file=output_prefix+'pickled_dmp.parsed_blast'; verbose('candidates_parsed_tblastx_file : '+candidates_parsed_tblastx_file+'\n')
  pretty_blast_candidates_file=output_prefix+'blastx.candidates.readable'; verbose('pretty_blast_candidates_file : '+pretty_blast_candidates_file+'\n')
  blasts_with_gene_in_the_middle=output_prefix+'blastx.gene_in_the_middle';verbose('blasts_with_gene_in_the_middle : '+blasts_with_gene_in_the_middle+'\n')
  exonerate_folder=output_prefix+'exonerate';                              verbose('exonerate_folder : '+exonerate_folder+'\n')
  dumped_output_filename=output_prefix+'pickled_dmp.outputs';                  verbose('dumped_output_filename : '+dumped_output_filename+'\n')
  gff_output_file=output_prefix+'output.gff';                              verbose('gff_output_file : '+gff_output_file+'\n')
  fasta_output_file=output_prefix+'output.fa';                             verbose('fasta_output_file : '+fasta_output_file+'\n')
  cds_output_file=output_prefix+'output.cds.fa';                           verbose('cds_output_file : '+cds_output_file+'\n')
  pretty_output_file=output_prefix+'output.pretty';                        verbose('pretty_output_file : '+pretty_output_file+'\n')
  secis_output_file =output_prefix+'output.secis.fa';                      verbose('secis_output_file : '+secis_output_file+'\n')
  all_secis_output_file =output_prefix+'output.all_secis';                 verbose('all_secis_output_file : '+all_secis_output_file+'\n')
  all_secis_output_fasta_file =output_prefix+'output.all_secis.fasta';     verbose('all_secis_output_fasta_file : '+all_secis_output_fasta_file+'\n')
  all_secis_output_gff_file= output_prefix+'output.all_secis.gff';         verbose('all_secis_output_gff_file : '+all_secis_output_gff_file+'\n')
  filtered_secis_output_file =output_prefix+'output.filtered_secis';       verbose('filtered_secis_output_file : '+filtered_secis_output_file+'\n')
  secis_images_folder   =output_prefix+'output.secis.images';              verbose('secis_images_folder  : '+secis_images_folder+'\n')

  #write( "******* " + str(load_chromosome_lengths), 1)

  ## printing header
  write(' '+'_'*54+' \n '+'_|*|'*4, how='yellow'); write('      Seblastian      '); write('|*|_'*4+' \n', 1, how='yellow')
  non_def_options_str= join(  [ '-'+str(k)+' '+str(opt[k])  for k in sorted(opt.keys()) if def_opt.has_key(k) and  def_opt[k] != opt[k]], '\n')
  write('      current non-default options:\n'+non_def_options_str+'\n', 1)
  # setting hierarchy: if one specify routine 1, also 2 and 3 must be executed 
  hierarchy_options=['search', 'load_secis', 'extract', 'find_uga', 'blast', 'fetch', 'parse', 'exonerate', 'filter', 'output'] 
  next_must_be_active=0
  for i in hierarchy_options: 
    if opt[i]: next_must_be_active=1
    elif next_must_be_active: opt[i]=1
  if opt['only_search']: opt['out_all_secis']=1
  #setting colored keyword. by now they are set as in selenoprofiles
  colored_keywords={'WARNING':'red', 'ERROR': 'red,underscore', '*OK':'green,reverse', '__TIME':'cyan'}
  for o in hierarchy_options: colored_keywords[upper(o)]='yellow'
  for l in output_labels: colored_keywords[l]='green'
  set_MMlib_var('colored_keywords', colored_keywords)
  # computing chromosome lengths if necessary, loading them
  chromosome_length_file=output_prefix+'chromosome_lengths'
  if not is_file(chromosome_length_file):
    service('Computing length of all chromosomes --> '+chromosome_length_file)
    bbash('fastalength '+target_file+ ' > '+temp_folder+'chrom_lengths && mv '+temp_folder+'chrom_lengths '+chromosome_length_file)
                            
  write('  Loading length of all chromosomes from: '+chromosome_length_file, 1)
  max_length=63
  if opt['fix']:   max_length=1000
  load_chromosome_lengths(chromosome_length_file,  max_length) # this will make the chromosome lengths available to selenoprofiles module and to MMlib; it's also checking the maximum number of characters for the fasta headers
  global chromosome_lengths; chromosome_lengths = get_MMlib_var('chromosome_lengths')
  
  if opt['time']:  #setting vars for timing execution
    import time
    global times_list; times_list=[]
    def time_report(name=''):
      """ compute and print a report of the time spent for an operation, computing it as the difference between this call of time.clock() and the last one. Uses variable times_list to keep track of times between calls.
      before the first printing time_report is called, it must be initialized by a call before the first operation is run. this won't print anything but it will fill the first value of times_list.
      if no name is specified, silent mode is used:  nothing will be printed and the current time.clock() will be appended to time_list """
      times_list.append(time.clock())
      if len(times_list)>1 and name:        write('__TIME for '+name+' : '+str( times_list[-1]-times_list[-2] )+' s ', 1)

  ## important global variables
  global secis_hash; secis_hash={} #index: secis  (index is a int)
  global blast_candidates_hash; blast_candidates_hash={} #index: SeB_blast
  global exonerate_candidates_hash; exonerate_candidates_hash={} #index: SeB_exonerate
  global output_candidates_hash; output_candidates_hash={} # label: {index: SeB_blast or SeB_exonerate}
  for label in output_labels: output_candidates_hash[label]={}

###########################################################
################ PROGRAM START
  verbose('# PROGRAM START! \n')
  try:
  #######################################################################################################################################  
    ### running SECIS search (infernal or others)
    all_secis=[]
    if opt['load_secis'] or not is_file(dumped_secis_filename):
      if 'infernal' in active_methods:
        if opt['search'] or not is_file(infernal_output_file):
          if opt['time']: time_report()
          secis_filter_cmnd='';       single_processor_infernal= (processors==1) or (opt['infernal_no_mpi'])
          if not opt['infernal_score'] is None:  secis_filter_cmnd+=' -T '+str(opt['infernal_score'])
          if not opt['infernal_evalue'] is None: secis_filter_cmnd+=' -E '+str(opt['infernal_evalue'])
          write( "SEARCH : Running infernal with model "+infernal_cm+' '+'('*int(bool(secis_filter_cmnd))+secis_filter_cmnd+' )'*int(bool(secis_filter_cmnd))+ ' --> '+infernal_output_file+'\n')
          b=[0, ''] #just to pass checks later
          if not single_processor_infernal:
            cmnd='mpirun -c '+str(processors)+' cmsearch --mpi '+secis_filter_cmnd+' --toponly '*(int(bool(opt['no_complement']))) +' '+infernal_cm+' '+target_file+' > '+temp_folder+'temp_infernal && mv '+temp_folder+'temp_infernal '+infernal_output_file;    b=bash(cmnd)
          #parsing errors if multithread is not enabled with infernal
          if b[0] and 'Failed to parse command line: No such option "--mpi".' in b[1]:
            write('WARNING infernal cmsearch could not be run with multiprocessor option. Please reinstall it using ./configure --enable-mpi', 1);       single_processor_infernal=True
          elif b[0] and 'not found' in b[1]:
            write('WARNING infernal cmsearch could not be run with multiprocessor option. mpirun was not found', 1);       single_processor_infernal=True
          elif b[0]: raise Exception, "check your infernal installation. You may want to try Seblastian option -infernal_no_mpi. ERROR with command *** "+cmnd+' *** :'+b[1]
          if single_processor_infernal: 
            bbash('cmsearch '+secis_filter_cmnd+' --toponly '*(int(bool(opt['no_complement']))) +' '+infernal_cm+' '+target_file+' > '+temp_folder+'temp_infernal && mv '+temp_folder+'temp_infernal '+infernal_output_file)
          if opt['time']: time_report('infernal search')
        ## loading SECIS from infernal
        write('  Loading SECIS from '+infernal_output_file+' ... ')  
        all_infernal_secis=parse_infernal(infernal_output_file).all()
        all_infernal_secis.sort(key=lambda x:x.score, reverse=True)  ## sorting secis by score
        for infernal_obj in all_infernal_secis:
          if (opt['infernal_score'] is None or infernal_obj.score>opt['infernal_score']) and (opt['infernal_evalue'] is None or infernal_obj.evalue<opt['infernal_evalue']): #filtering score/evalue
            infernal_obj.cmali= infernal_model_alignment
            infernal_obj.target= target_file
            if opt['n_secis'] and len(all_secis) >= opt['n_secis']: break
            x=SeB_secis(infernal_obj)
            all_secis.append(  x  )
        write('number: '+str(len(all_secis)), 1)
        del all_infernal_secis

      if 'secisearch1' in active_methods:
        if opt['search'] or not is_file(secisearch_output_file):
          if opt['time']: time_report()
          write( "SEARCH : Running secisearch with pattern "+secisearch_pattern+' --> '+secisearch_output_file+'\n')
          bbash('SECIScript.pl -f '+target_file+'  -p '+secisearch_pattern+'  -complement '*(int(not bool(opt['no_complement'])))+' -temp '+temp_folder+'  >  '+temp_folder+'temp_secisearch && mv '+temp_folder+'temp_secisearch '+secisearch_output_file)
          if opt['time']: time_report('secisearch1')      
        ## loading SECIS from secisearch
        write('  Loading SECIS from '+secisearch_output_file+' ... ')
        cline= bbash('head -1 '+secisearch_output_file)
        all_secisearch_secis=[]
        if not 'Sorry, no SECIS elements were found! \n' in  cline:                                                            
          for secisearch_obj in parse_secisearch(secisearch_output_file):                     all_secisearch_secis.append( secisearch_obj )
          write(str(len(all_secisearch_secis)), 1)
          for secisearch_obj in all_secisearch_secis:          
            secisearch_obj.target= target_file
            x=SeB_secis(secisearch_obj)            
            all_secis.append( x )
        else: write('0', 1)
        del all_secisearch_secis
      
      if 'covels' in active_methods:
        if opt['time']: time_report()
        if opt['search'] or not is_file(covels_output_file):
          secis_filter_cmnd=''; 
          if not opt['covels_score'] is None: secis_filter_cmnd+=' -t '+str(opt['covels_score']) #can't be None actually!!
          write( "SEARCH : Running covels with model "+covels_cm+' '+'('*int(bool(secis_filter_cmnd))+secis_filter_cmnd+' )'*int(bool(secis_filter_cmnd))+' --> '+covels_output_file+'\n')
          bbash('covels-SE '+secis_filter_cmnd+' -c '*(int(not bool(opt['no_complement'])))+'  '+covels_cm+' '+target_file+' >  '+temp_folder+'temp_covels && mv '+temp_folder+'temp_covels '+covels_output_file )
          if opt['time']: time_report('covels search')      
        write('  Loading SECIS from '+covels_output_file+' ... ')
        all_covels_secis=[ covels_obj     for covels_obj in parse_covels(covels_output_file)  if covels_obj.score >= opt['covels_score'] ] #to allow using covels score even when loading.
        write(str(len(all_covels_secis)), 1)
        for covels_obj in all_covels_secis:          
          covels_obj.target= target_file
          x=SeB_secis(covels_obj)
          all_secis.append( x )
        del all_covels_secis

      """ created for testing purpose. You can load SECIS from a GFF as well"""
      if 'gff' in active_methods:
        if not is_file(secis_gff_infile): raise Exception, "ERROR file "+secis_gff_infile+" could not be found! "
        write('  Loading SECIS from '+secis_gff_infile+' ... ')
        all_secis_from_gff = load_all_genes(secis_gff_infile, uniqid_function=lambda x:x.split('\t')[-1].split()[0], tag='*')
        write(str(len(all_secis_from_gff)), 1)
        for gene_obj in all_secis_from_gff:      
          gene_obj.target=target_file
          all_secis.append( SeB_secis(gene_obj) )
        del all_secis_from_gff

      def choose_secis(secis1, secis2):
        """ method to keep best secis when 2 overlaps. it chooses depending on the source of prediction. it also keep track of all predictions merged keeping the best scores (or pattern) for ecah type"""
        if secis2.infernal_score > secis1.infernal_score:
          secis2.covels_score= max(  [secis1.covels_score, secis2.covels_score]);    secis2.secisearch_pattern=max([secis1.secisearch_pattern, secis2.secisearch_pattern])
          return secis2
        elif secis2.infernal_score < secis1.infernal_score:
          secis1.covels_score= max(  [secis1.covels_score, secis2.covels_score]);    secis1.secisearch_pattern=max([secis1.secisearch_pattern, secis2.secisearch_pattern])
          return secis1

        elif secis2.covels_score > secis1.covels_score:
          secis2.infernal_score= max(  [secis1.infernal_score, secis2.infernal_score]);    secis2.secisearch_pattern=max([secis1.secisearch_pattern, secis2.secisearch_pattern])
          return secis2
        elif secis2.covels_score < secis1.covels_score:
          secis1.infernal_score= max(  [secis1.infernal_score, secis2.infernal_score]);    secis1.secisearch_pattern=max([secis1.secisearch_pattern, secis2.secisearch_pattern])
          return secis1

        elif secis2.secisearch_pattern > secis1.secisearch_pattern:
          secis2.infernal_score= max(  [secis1.infernal_score, secis2.infernal_score]);    secis2.covels_score= max(  [secis1.covels_score, secis2.covels_score])
          return secis2
        elif secis2.secisearch_pattern < secis1.secisearch_pattern:
          secis1.infernal_score= max(  [secis1.infernal_score, secis2.infernal_score]);    secis1.covels_score= max(  [secis1.covels_score, secis2.covels_score])
          return secis1
        return secis1
        
      #merging, removing the same secis predicted by different methods
      if len(active_methods)>1:
        write('  Merging SECIS from different methods... ')
        all_secis =    merge_genes(all_secis, phase=False, inplace=False, mode=choose_secis)   
        write('number after merging: '+str(len(all_secis)), 1)
      ### filtering SECIS elements!
      if not opt['no_secis_filter']: 
        if opt['out_filtered_secis']:
          filtered_secis_output_file_h=open(temp_folder+'out_filtered_secis', 'w')
        write('  Filtering SECIS: '+str(len(all_secis))+'  -->  ')
        out_s=[]
        for x in all_secis:
          a=secis_filter(x) #filtering secis
          if a[0]: out_s.append(x)
          elif opt['out_filtered_secis']:
            print >> filtered_secis_output_file_h, x.summary(add=' Filtered out because: '+a[1])
        if opt['out_filtered_secis']: filtered_secis_output_file_h.close()
        all_secis=out_s
        write(len(all_secis), 1)
        if opt['out_filtered_secis']:
          bbash('mv '+temp_folder+'out_filtered_secis '+filtered_secis_output_file)
          write('  Wrote filtered out SECIS elements in: '+filtered_secis_output_file, 1)

        
      #excluding SECIS in -exclude
      if opt['exclude']:
        write('  Excluding SECIS found in file: '+opt['exclude']+' ')
        all_secis_to_exclude=[]
        for exclude_file in exclude_files:
          extension=exclude_file.split('.')[-1]
          if      extension == 'gff':      all_secis_to_exclude.extend(  load_all_genes(exclude_file, tag='*', uniqid_function=lambda x:x.split('\t')[-1].split()[0])   )
          elif    extension == 'infernal':     all_secis_to_exclude.extend(  parse_infernal(exclude_file).all()  )
          elif    extension == 'covels':       all_secis_to_exclude.extend(  parse_covels(exclude_file).all()  )
          elif    extension == 'secisearch1':   all_secis_to_exclude.extend(  parse_secisearch(exclude_file).all()  )
          elif    extension == 'erpin':        all_secis_to_exclude.extend(  parse_erpin(exclude_file).all()  )
          else:   raise Exception, "ERROR file type not recognized for: "+exclude_file
        for s in all_secis_to_exclude:      s.id='exclude' #tagging the secis to remove them later
        write('('+str(len(all_secis_to_exclude))+') ... ')
        def keep_secis_in_excluded_list(secis1, secis2):
          """ sounds counterintuitive, but this is the function fed to merge_genes to keep always these genes when checking overlaps. so then we can just remove all these gene since we know their pecularity, the id"""
          if    secis1.id=='exclude': return secis1
          elif  secis2.id=='exclude': return secis2
          else: raise Exception, "ERROR excluding! "
        all_secis =    merge_genes(all_secis, phase=False, inplace=False, mode=keep_secis_in_excluded_list)
        for i in range( len(all_secis)-1, -1, -1 ):
          if all_secis[i].id=='exclude': all_secis.pop(i)
        write('secis number after excluding: '+str(len(all_secis)), 1)

      # filling secis_hash
      for secis_index, secis_obj in enumerate ( all_secis ): #here secis_index is 0 based, but it will 1 based in the next assignments
        secis_obj.target=target_file
        secis_obj.id=str(secis_index+1)
        secis_hash[secis_index+1]=secis_obj
      dumped_secis_filename_h=open(dumped_secis_filename, 'w')
      write('  Saving processed SECIS elements --> '+dumped_secis_filename, 1)
      pickle.dump(secis_hash, dumped_secis_filename_h)
      dumped_secis_filename_h.close()
    else: 
      write('  Loading all SECIS from saved file: '+dumped_secis_filename)
      secis_hash=pickle.load( open(dumped_secis_filename, 'r'))
      write(' (tot:'+str(len(secis_hash))+')', 1)
      if opt['fix']:

        for i in secis_hash:
          s=secis_hash[i]
          s.refined_ss= replace(s.refined_ss, '_', '-')

        dumped_secis_filename_h=open(dumped_secis_filename, 'w')
        write('  Saving processed SECIS elements --> '+dumped_secis_filename, 1)
        pickle.dump(secis_hash, dumped_secis_filename_h)
        dumped_secis_filename_h.close()

      
    write_all_secis=False; write_all_secis_fasta=False
    if opt['out_all_secis']:
      write_all_secis=        not is_file(all_secis_output_file) or opt['output']
      write_all_secis_fasta=  opt['out_secis_fasta'] and (not is_file(all_secis_output_fasta_file) or opt['output'])
      write_all_secis_gff  =  opt['out_gff'] and   (not is_file(all_secis_output_gff_file) or opt['output'])
      if opt['out_secis_image']:       secis_images_folder=Folder(secis_images_folder) #creating folder is necessary      
      if write_all_secis or write_all_secis_fasta or write_all_secis_gff or opt['out_secis_image']:
        if write_all_secis:
          all_secis_output_file_h=open(all_secis_output_file, 'w')
          write('  Writing all SECIS to file --> '+all_secis_output_file, 1)
        if write_all_secis_fasta:
          all_secis_output_fasta_file_h=open(all_secis_output_fasta_file, 'w')
          write('  Writing all SECIS in fasta to file --> '+all_secis_output_fasta_file, 1)
        if write_all_secis_gff:
          all_secis_output_gff_file_h=open(all_secis_output_gff_file, 'w')
          write('  Writing all SECIS in gff to file --> '+all_secis_output_gff_file, 1)
        if opt['out_secis_image']:
          write('  Writing all SECIS images in folder: '+secis_images_folder, 1)

        for secis_index in range(1, len(secis_hash)+1):
          secis_obj=secis_hash[secis_index]
          if write_all_secis:          print >> all_secis_output_file_h, secis_obj.summary()+'\n'
          if write_all_secis_fasta:    print >> all_secis_output_fasta_file_h, secis_obj.fasta()
          if write_all_secis_gff:      print >> all_secis_output_gff_file_h, secis_obj.gff()
          secis_image_filename=         secis_images_folder+'SECIS.'+str(secis_index)+'.png'
          secis_image_refined_filename= secis_images_folder+'SECIS.'+str(secis_index)+'.refined.png'
          if opt['out_secis_image'] and (not is_file(secis_image_filename) or opt['output']):            secis_obj.draw_image(secis_image_filename)

        if write_all_secis:         all_secis_output_file_h.close() 
        if write_all_secis_fasta:   all_secis_output_fasta_file_h.close()
        if write_all_secis_gff:     all_secis_output_gff_file_h.close()

    if opt['only_search']: 
      write('Finished', 1)
      raise nohitsException
    if not secis_hash: 
      write('No SECIS found, exiting...', 1)
      raise nohitsException

  #######################################################################################################################################
   ### extracting upstream sequences
    if opt['extract'] or not is_file(upstream_seqs_file):
      if opt['time']: time_report()
      write( "EXTRACT : Obtaining the sequences upstream of all SECISes --> "+upstream_seqs_file, 1 )
      wrote_at_least_one_upstream_seq=False
      for secis_index in range(1, len(secis_hash)+1):
        secis_obj=secis_hash[secis_index]      
        verbose("extracting seq upstream of SECIS."+secis_obj.id+'\n')
        length_upstream_used_for_this_secis = length_upstream #it may change if upstream region is out of bounds of the chromosome
        upstream_gene_obj=secis_obj.upstream(length_gap, length_upstream)
        if max(upstream_gene_obj[0])>1: 
          chrom_file = fastafetch(split_folder, secis_obj.chromosome, target_file)
          if secis_obj.strand=='+':
            if upstream_gene_obj[0][0]<=0: #out of bounds on the left
              length_upstream_used_for_this_secis=length_upstream + upstream_gene_obj[0][0] -1
              upstream_gene_obj.exons[0][0]=1
            fastasubseq(chrom_file, upstream_gene_obj[0][0]-1, length_upstream_used_for_this_secis, temp_folder+'fragment_upstream')
          elif secis_obj.strand=='-':
            try:
              length_upstream_used_for_this_secis = fastasubseq(chrom_file, upstream_gene_obj[0][0]-1, length_upstream, temp_folder+'to_complement')
              bbash('fastarevcomp '+temp_folder+'to_complement > '+temp_folder+'fragment_upstream')
            except Exception: #out of bounds
              length_upstream_used_for_this_secis=0           
          if length_upstream_used_for_this_secis >= 15:
                ## here writing sequence, it passed the "out of bounds" tests
            wrote_at_least_one_upstream_seq=True
            bbash('echo ">upstream_of_SECIS_'+secis_obj.id+' length='+str(length_upstream_used_for_this_secis)+'" >> '+temp_folder+'temp_upstream_seqs')
            bbash("gawk 'NR>1' " + temp_folder+'fragment_upstream >> '+temp_folder+'temp_upstream_seqs')
        else:            verbose('SECIS.'+secis_obj.id+' out of boundaries\n')
      #moving file to final destination. so it computation is interrupted, no partial file is produced
      if not wrote_at_least_one_upstream_seq: 
        write('No SECIS had regular upstream sequence, exiting...', 1)
        raise nohitsException
      bbash('mv  '+temp_folder+'temp_upstream_seqs '+upstream_seqs_file)
      if opt['time']: time_report('extract phase')

  #######################################################################################################################################
   ### writing uga positions in uga_positions file. format:   upstream_of_SECIS_n seq_length pos1 pos2 pos3
    if opt['find_uga'] or not is_file(uga_positions_file):
      if opt['time']: time_report()
      write( "FIND_UGA : Finding all UGA codons in the upstream sequences --> "+uga_positions_file +'\n')
      uga_positions_file_h = open(temp_folder+'temp_uga_pos', 'w')
      for title, seq in parse_fasta(upstream_seqs_file):
        if 1:
          seq=upper(seq)
          print >> uga_positions_file_h, title, 
          f=seq.find('TGA', 0)
          while f!=-1:
            print >> uga_positions_file_h, f+1,
            f=seq.find('TGA', f+1)
          print >> uga_positions_file_h, ""
      uga_positions_file_h.close()
      #moving file to final destination. so it computation is interrupted, no partial file is produced
      bbash('mv  '+temp_folder+'temp_uga_pos '+uga_positions_file)
      if opt['time']: time_report('find_uga phase')

  ######################################################################################################################################
   ### running blastx and formatting output
    if opt['blast'] or not is_file(blastx_output_file):
      if opt['time']: time_report()
      if opt['net']:
        #running netblast #### having problem with this! doens't work!
        write( "BLAST : Running netblastx at ncbi searching against db "+prot_database+ " --> "+blastx_output_file+'\n' )
        bbash('blastcl3 -b 5000 -p blastx -S 1 -I  -i '+upstream_seqs_file+' -d '+prot_database+' -e '+blast_evalue+' '+blast_opt+' > '+temp_folder+'temp_blastx && mv '+temp_folder+'temp_blastx '+blastx_output_file)
      else:
        #running local blastall instead (default)
        write( "BLAST : Running blastx search against db "+prot_database+ " --> "+blastx_output_file+'\n' )
        bbash('blastall -b 5000 -p blastx -S 1 -I  -i '+upstream_seqs_file+' -d '+prot_database+' -a '+str(processors)+' -e '+blast_evalue+' '+blast_opt+' > '+temp_folder+'temp_blastx && mv '+temp_folder+'temp_blastx '+blastx_output_file)
      if opt['time']: time_report('blast phase')
    if opt['blast'] or not is_file(parsed_blastx_output_file):
      if opt['time']: time_report()
      write( "  Formatting blastx output with blaster_parser.g --> "+parsed_blastx_output_file+'\n' )
      bbash('blaster_parser.g -v TLENGTH=1 -v QLENGTH=1 '+blastx_output_file+' > '+temp_folder+'temp_blastx.parsed && mv '+temp_folder+'temp_blastx.parsed '+parsed_blastx_output_file)
      if opt['time']: time_report('blaster parsing')
      if not bbash('head '+parsed_blastx_output_file): 
        write('No blast hits found, exiting...', 1)
        raise nohitsException

  #######################################################################################################################################
   ### annotating proteins hit which contains a U
    if opt['fetch'] or not is_file(fetched_sequences_gis_file):
      if opt['time']: time_report()
      write( "FETCH : Fetching all protein sequences that hit upstream of any secis, keeping those with at least a U  --> " +fetched_sequences_gis_file +'\n')
      bbash("gawk '{print $1}' "+parsed_blastx_output_file+' | sort | uniq > ' +temp_folder+'all_target_protein_gis')    
      if not opt['net']:
        bbash('fastacmd -d '+prot_database+' -i '+temp_folder+'all_target_protein_gis'+' > '+temp_folder+'all_target_protein_gis.fasta')
      else:
        fetched_hash = fetch_ncbi_sequences.main({'silent':1, 'f':1, 'm':'P', 'i':temp_folder+'all_target_protein_gis'})
        writing_file_h=open(temp_folder+'all_target_protein_gis.fasta', 'w')
        for gi_id in fetched_hash:        print >> writing_file_h, ">"+fetched_hash[gi_id][0]+'\n'+fetched_hash[gi_id][1]
        writing_file_h.close()
      bbash('gawk \'{  if ( /^>/) { if (print_this) { printf "%s\\n", substr(title, 2) };	title=$1; print_this=0 }  else { if (/U/){print_this=1}  } }END{if (print_this){printf "%s\\n", substr(title, 2)}}\' ' +temp_folder+'all_target_protein_gis.fasta  > '+temp_folder+'temp_fetched_sequences_gis_file && mv '+temp_folder+'temp_fetched_sequences_gis_file '+fetched_sequences_gis_file  )
      if opt['time']: time_report('fetch phase')

  #######################################################################################################################################
   ### parsing blast outputs, filling blast_candidates_hash (and blast_candidates_indexes with their ids converted to int)
    blast_candidates_indexes=[]
    category_hierarchy=['ALI_U','QUE_U', 'ALI_X', 'ALI_C', 'END']
    # loading in memory the ids of proteins with U
    global proteins_with_annotated_u; proteins_with_annotated_u={}
    for l in open(fetched_sequences_gis_file, 'r'):      proteins_with_annotated_u[l[:-1]]=1
    if opt['parse'] or not is_file(candidates_parsed_tblastx_file):
      if opt['time']: time_report()
      summary_hash={}
      parsed_blastx_output_file_h=open(parsed_blastx_output_file, 'r')
      write( "PARSE : Looking for UGA to U/C/X alignments or hits to known selenoproteins or complete CDS upstreams of SECIS --> "+candidates_parsed_tblastx_file+'\n' )
      pretty_blast_candidates_file_h=open(temp_folder+'temp_pretty_blast_candidates_file', 'w')
      if opt['check_middle']:        blasts_with_gene_in_the_middle_h=open(temp_folder+'blasts_with_gene_in_the_middle', 'w')
      #reading and parsing the first blast line
      bline=parsed_blastx_output_file_h.readline()
      if bline:       blasth=SeB_blast(bline)

      for secis_index in  range(1, len(secis_hash)+1):       
        #service('parsing blastx for secis '+str(secis_index)+' ...')
        #cycling each secis... 
        if not bline: break
        while bline and int(blasth.id) < secis_index: #skipping the useless ones, getting into the first blast hit of this secis
          bline=parsed_blastx_output_file_h.readline()
          if bline:       blasth=SeB_blast(bline)
        if int(blasth.id) > secis_index: continue #skipping secis index for which no blast hit at all is present
        ## loading uga positions
        c=bbash("gawk '$1==\"upstream_of_SECIS_"+str(secis_index)+"\"' "+uga_positions_file).split()
        uga_codons_hash={} #will keep track of the positions 1 based of the G of TGA (choosing central pos so I don't have to split + and - )
        if len(c)> 2 : #at least one UGA position
          for i in c[2:]:
            uga_codon=blasth.blast_query.subseq(int(i), 3)
            uga_codons_hash[ uga_codon[0][0]+1 ]=True
        verbose('SECIS: '+str(secis_index)+' ; upstream UGA positions: '+str( len(uga_codons_hash) )+'\n' )
        all_blast_hit_for_this_secis=[]
        while bline and  int(blasth.id) == secis_index:
          ## parsing each blast hit (those corresponding to the secis element we're parsing)
          ### determining its category, that is to say whether we see a U aligned to UGA, or a C etc
          x=blasth
          all_blast_hit_for_this_secis.append(blasth)
          if opt['blast_filter'] is None or eval(opt['blast_filter']): #applying filter to blast hits if active
            blasth.set_ugas(uga_codons_hash);   ### blasth.target=target_file
            blasth.assign_category()
            #saving blasth if it's the best so far for this secis: category is the major determinant. then, best evalues are parsed first, so we simply don't overwrite a blasth of equal category
            if blasth.category!='NOT' and (  not blast_candidates_hash.has_key(secis_index) or  category_hierarchy.index( blasth.category ) <  category_hierarchy.index( blast_candidates_hash[secis_index].category )  ):
              blast_candidates_hash[secis_index]=blasth
          bline=parsed_blastx_output_file_h.readline()
          if bline:  ##### added just to avoid a bug of blast. in a single case among all my runs, I have blast reporting a hit on the negative strand, which wasn't requested! 
            blasth=SeB_blast(bline)
            while blasth is None and bline: 
              bline=parsed_blastx_output_file_h.readline()
              blasth=SeB_blast(bline)
        
        ## outputing of this phase; format: abstract_line  for candidates_parsed_tblastx_file, summary for pretty
        if blast_candidates_hash.has_key(secis_index):  
          this_category=blast_candidates_hash[secis_index].category
          blast_hit_same_query_as_candidate=[]
          blast_candidates_hash[secis_index].is_candidate=True #adding property to recognize it later
          for b in all_blast_hit_for_this_secis:
            if  b.query.chromosome == blast_candidates_hash[secis_index].query.chromosome:            blast_hit_same_query_as_candidate.append(b) #including the candidate

          if opt['check_middle']:          all_blast_hit_for_this_secis = merge_by_colinearity_p2g_with_different_queries( all_blast_hit_for_this_secis ) #also running set_ugas on the hits joined with our candidates
          else:                            all_blast_hit_for_this_secis = merge_by_colinearity_p2g_with_different_queries( blast_hit_same_query_as_candidate ) #NB replacing all_blast_hit_for_this_secis with joined blast_hit_same_query_as_candidate, HAVE TO OPTIMIZE IF WE DROP check_middle option

          for b in all_blast_hit_for_this_secis:          
            if hasattr(b, 'is_candidate'):
              blast_candidates_hash[secis_index]=b #our candidate blast hit, maybe enriched with other exons
              blast_candidates_hash[secis_index].category=this_category
          genes_in_the_middle=[]
          if opt['check_middle']:
            ## make a gene containing the projection of our gene plus all genes that overlaps with it
            projection_all_possible_exons=blast_candidates_hash[secis_index].copy()
            for b in all_blast_hit_for_this_secis:
              if b!=blast_candidates_hash[secis_index] and blast_candidates_hash[secis_index].overlaps_with(b):
                projection_all_possible_exons=projection_all_possible_exons.union_with(b)
            ## check for genes in the middle:   any blast hit meeting a specific filter (option), downstream of our gene, and not overlapping the gene above        
            for b in all_blast_hit_for_this_secis:
              if b.is_downstream_of( blast_candidates_hash[secis_index] ) and not b.overlaps_with( projection_all_possible_exons  ) and genes_in_the_middle_filter(b, blast_candidates_hash[secis_index]): 
                genes_in_the_middle.append(b)

          if genes_in_the_middle: # there are one or more genes between the blast hit and the secis
            print >> blasts_with_gene_in_the_middle_h, "*"*150+'\n'+str(blast_candidates_hash[secis_index])+'\n'+'genes in the middle'.center(150, '^')+join([str(a) for a in genes_in_the_middle], '\n')
          else:
            ### blast hit accepted!
            blast_candidates_hash[secis_index].remove_Xs()
            blast_candidates_hash[secis_index].fill_full_query_information() # to be sure we don't have to recompute it when we load the blast it from dumped results
            print >> pretty_blast_candidates_file_h, blast_candidates_hash[secis_index].summary() 
            blast_candidates_indexes.append(secis_index)
            summary_hash.setdefault(blast_candidates_hash[secis_index].category, 0); summary_hash[blast_candidates_hash[secis_index].category]+=1

      pretty_blast_candidates_file_h.close()
      if opt['check_middle']:    blasts_with_gene_in_the_middle_h.close()
      parsed_blastx_output_file_h.close()

      write('  Saving parsed blastx candidates --> '+candidates_parsed_tblastx_file, 1 )
      candidates_parsed_tblastx_file_h=open(candidates_parsed_tblastx_file, 'w')
      pickle.dump( blast_candidates_hash,  candidates_parsed_tblastx_file_h )
      candidates_parsed_tblastx_file_h.close()

      bbash('mv '+temp_folder+'temp_pretty_blast_candidates_file '+pretty_blast_candidates_file)
      if opt['check_middle']:        bbash('mv '+temp_folder+'blasts_with_gene_in_the_middle '+blasts_with_gene_in_the_middle)
      summary_categories=join([ k+': '+str(summary_hash[k]) for k in sorted(summary_hash.keys(), key=category_hierarchy.index) ], ' ')
      write('  You can check the candidates in human-readable format in file: '+pretty_blast_candidates_file+' ( '+summary_categories+' ) ', 1)
      if opt['time']: time_report('parse phase')
    else:
      write('  Loading blast candidates from saved file: '+candidates_parsed_tblastx_file+' ') # from pickled file
      candidates_parsed_tblastx_file_h=open(candidates_parsed_tblastx_file, 'r')
      blast_candidates_hash = pickle.load( candidates_parsed_tblastx_file_h )
      candidates_parsed_tblastx_file_h.close()
      summary_hash={}    
      for secis_index  in sorted(  blast_candidates_hash.keys() ):
        blasth = blast_candidates_hash[secis_index]
        blast_candidates_indexes.append(secis_index)
        summary_hash.setdefault(blasth.category, 0); summary_hash[blasth.category]+=1
        blasth.secis=get_secis(secis_index)
        if blasth.chromosome != blasth.secis.chromosome: raise Exception, "ERROR the blast candidate and its secis are not on the same chromosome!"
        
        
      summary_categories=join([ k+':'+str(summary_hash[k]) for k in sorted(summary_hash.keys(), key=category_hierarchy.index) ], ' ')+' (tot:'+str(len(blast_candidates_hash))+')'
      write( summary_categories , 1)

  #######################################################################################################################################
    ### Exonerate phase: obtaining complete gene predictions
    def run_cyclic_exonerate(blasth, outfile):
      """ Run the cyclic exonerate routine imported from selenoprofiles and returns an exonerate object """
      #fetching the whole seq of the original blast target (protein)
      exonerate_query_full_title = blasth.get_query_full_title()
      exonerate_query_seq = blasth.get_query_full_sequence()
      exonerate_query_seq = replace(exonerate_query_seq, 'J', 'X')
      if exonerate_query_full_title.startswith('lcl|'): exonerate_query_full_title=exonerate_query_full_title[4:]
      query_profile= alignment() # fake alignment, it's actually just one sequence
      query_profile.add(exonerate_query_full_title, exonerate_query_seq)
      #print query_profile
      
      exonerate_options={}
      if opt['genetic_code']:        exonerate_options['-geneticcode']=str(opt['genetic_code'])

      
      spexoneratehit=cyclic_exonerate(query_profile, target_file, outfile=outfile, seed=blasth, extension=cyclic_exonerate_ext, exonerate_options=exonerate_options)
      return SeB_exonerate( exonerateh=spexoneratehit )     #running cyclic_exonerate as imported from selenoprofiles, converting to SeB_classc

    def keep_blast(exonerateh):
      """ Decides if the exonerateh should be discarded in favor of the blast seed"""
      if exonerateh.label in keep_blast_labels: return True
      elif exonerateh.category in ['ALI_U', 'ALI_C', 'ALI_X'] and exonerateh.label in ['no_sec_candidate']: return True
      return False
      
    if not (opt['exonerate'] or opt['filter'] ) and is_file(dumped_output_filename): 
      write('  Loading all filtered results from saved file: '+dumped_output_filename, 1)
      dumped_output_filename_h=open(dumped_output_filename)
      output_candidates_hash = pickle.load(dumped_output_filename_h)
      dumped_output_filename_h.close()
      for label in output_candidates_hash:
        for secis_index in output_candidates_hash[label]:  output_candidates_hash[label][secis_index].secis=get_secis(secis_index)
    else:
      #parsing candidates running or loading each exonerate 
      if opt['time']: time_report()
      exonerate_folder=Folder(exonerate_folder) #creating folder if necessary, adding "/" if necessary
      write( "EXONERATE :  completing selenoprotein gene structures  --> "+exonerate_folder +'\n')    
      for secis_index in blast_candidates_indexes:
        blasth=blast_candidates_hash[secis_index]
        blasth.alignment.set_sequence('q',   replace(blasth.alignment.seq_of('q'), 'J', 'X') )
        exonerate_output_file=exonerate_folder+'SECIS.'+blasth.id+'.exout'
        write('-ID:'+str(secis_index).ljust(10)+(' blastx:'+blasth.category).ljust(15))
        if opt['exonerate'] or not is_file(exonerate_output_file):
          #running
          write(' -R> ')
          opt['filter']=1; opt['output']=1
          exonerateh= run_cyclic_exonerate(blasth, exonerate_output_file)
        else: 
          #loading
          write(' -L> ')
          spexoneratehit = superexoneratehit(); spexoneratehit.load(exonerate_output_file, seed=blasth)
          exonerateh     = SeB_exonerate(exonerateh=spexoneratehit)
        exonerateh.id=str(secis_index);   exonerateh.category=blasth.category;   exonerateh.secis=get_secis(secis_index); exonerateh.target=target_file
        exonerateh.query_full_seq  = blasth.query_full_seq; exonerateh.query_full_title = blasth.query_full_title; exonerateh.seed=blasth
        exonerateh.set_ugas(); 
        exonerateh.query_prot_length=blasth.query_prot_length
        exonerateh.assign_label()
        exonerate_candidates_hash[secis_index]=exonerateh
        write( ('exonerate:'+exonerateh.label).ljust(38) )        

        #if exonerate is empty or not consistent with blast, let's keep the information in the blast instead.
        if keep_blast( exonerateh ):        
          blasth.assign_label()
          outp2g=blasth         
          write( (' |->back2blastx:'+blasth.label).ljust(42)+' ')
        else:
          outp2g=exonerateh
          write(' '*43)

        ### filtering results
        #if opt['filter']:   ##it's True for sure
        x=outp2g
        if outp2g.label in output_labels:
          if not opt['custom_filter'] is None:         passed_filtering= eval(opt['custom_filter'])
          else:                                 passed_filtering= default_filter(x)
          if passed_filtering:
            output_candidates_hash[outp2g.label][secis_index]=outp2g
            write('*OK')
          else:  write('_no')
        else:                write('   ')
        write('|', 1)

      if opt['time']: time_report('exonerate/filtering phase (run or loaded)')
      #dumping output hash
      write('  Saving filtered results --> '+dumped_output_filename, 1) 
      dumped_output_filename_h=open(dumped_output_filename, 'w')
      pickle.dump(output_candidates_hash, dumped_output_filename_h)
      dumped_output_filename_h.close()

    ## output phase    
    #opening filehandlers
    if opt['time']: time_report()
    write_gff_file=opt['out_gff'] and (not is_file(gff_output_file) or opt['output'] )
    if write_gff_file:      gff_output_file_h=open(gff_output_file, 'w')
    write_fasta_file=opt['out_fasta'] and (not is_file(fasta_output_file) or opt['output'] )
    if write_fasta_file:    fasta_output_file_h=open(fasta_output_file, 'w')
    write_cds_file=opt['out_cds'] and (not is_file(cds_output_file) or opt['output'] )
    if write_cds_file:    cds_output_file_h=open(cds_output_file, 'w')
    write_pretty_file= opt['out_pretty'] and (not is_file(pretty_output_file) or opt['output'])
    if write_pretty_file:    pretty_output_file_h=open(pretty_output_file, 'w')
    write_secis_file=  opt['out_secis_fasta'] and (not is_file(secis_output_file) or opt['output'])
    if write_secis_file:     secis_output_file_h=open(secis_output_file, 'w')
  
    wrote_at_least_one_image=False
    if opt['out_secis_image']:       secis_images_folder=Folder(secis_images_folder) #creating folder is necessary      

    output_summary='( '+join(  [label+': '+str(len(output_candidates_hash[label]))+' ' for label in output_candidates_hash if output_candidates_hash[label]], '')+')'
    if output_summary=='( )': 
      write('No hits to output!', 1)
      if write_fasta_file:      fasta_output_file_h.close()
      if write_cds_file:        cds_output_file_h.close()
      if write_gff_file:        gff_output_file_h.close()
      if write_pretty_file:     pretty_output_file_h.close()
      if write_secis_file:      secis_output_file_h.close()
      raise nohitsException
    write("OUTPUT: printing information on results passing filtering "+output_summary, 1)
    def sorting_results_function(a, b):
      if a.secis.grade()<b.secis.grade():   return -1
      elif a.secis.grade()>b.secis.grade(): return 1
      else:         return cmp(int(a.id), int(b.id))
    for label in output_labels:
      p2g_outlist= sorted( output_candidates_hash[label].values(), cmp=sorting_results_function)
      for outp2g in p2g_outlist:
        secis_index=int(outp2g.id)
        outp2g=output_candidates_hash[label][secis_index]
        outp2g.secis.remove_Xs()
        write((outp2g.one_line_summary()).ljust(150)+' ' )
        if write_pretty_file: 
          write('>P ')
          print >> pretty_output_file_h, str(outp2g)
        if write_gff_file:    
          write('>G ')
          print >> gff_output_file_h,    outp2g.output_gff()
        if write_fasta_file:  
          write('>F ')
          print >> fasta_output_file_h,  outp2g.output_fasta()
        if write_cds_file:  
          write('>C ')
          print >> cds_output_file_h,  outp2g.output_cds()
        if write_secis_file:
          write('>S ')
          print >> secis_output_file_h, outp2g.output_secis()

        if opt['out_secis_image']:
          secis_image_filename= secis_images_folder+'SECIS.'+str(secis_index)+'.png'
          secis_image_refined_filename= secis_images_folder+'SECIS.'+str(secis_index)+'.refined.png'
          if not is_file(secis_image_filename) or opt['output']:   
            wrote_at_least_one_image=True
            outp2g.secis.draw_image(secis_image_filename)
            write('>I ')

        write('', 1)

    write('\n  Computation finished. These output files were produced:', 1)
    if write_fasta_file:      
      write(' '+fasta_output_file, 1)
      fasta_output_file_h.close()
    if write_cds_file:      
      write(' '+cds_output_file, 1)
      cds_output_file_h.close()
    if write_gff_file:        
      write(' '+gff_output_file, 1)
      gff_output_file_h.close()
    if write_pretty_file:     
      write(' '+pretty_output_file, 1)
      pretty_output_file_h.close()
    if write_secis_file:      
      write(' '+secis_output_file, 1)
      secis_output_file_h.close()
    if write_all_secis:           write(' '+all_secis_output_file, 1)
    if write_all_secis_fasta:     write(' '+all_secis_output_fasta_file, 1)
    if wrote_at_least_one_image:  write(' .png files in: '+secis_images_folder, 1)

    if opt['time']:       time_report('output phase')
    
    return  secis_hash, output_candidates_hash
  except nohitsException:    
    return  secis_hash, output_candidates_hash
    
###################
#### finished Main. here are some functions and the classes used :

def get_secis(secis_index):
    """ silly function"""
    if not secis_hash.has_key(secis_index):        raise Exception, "SECIS "+str(secis_index)+" not loaded!"
    return secis_hash[secis_index]

def secis_filter(self):
  """ Filter developed with Alexei to filter out unlikely structures. Returns (Result, description) where result is a boolean which has True for "good secis" and False for bad, and description is a string indicating the reason why it was filtered out ("ok" if it passes filtering)"""
  if not self.core or len(self.core[0])<2: return False, "no core aligned"
  core_seq_5 = join( [self.sequence()[ i-1  ] for i in self.core[0] ] , '' )
  core_seq_3 = join( [self.sequence()[ j-1  ] for j in self.core[1] ] , '' )
  if not 'GA' in core_seq_5 or not 'GA' in core_seq_3:   return False, "no GA-GA in core"

  a_length= self.apical_loop_length()
#  if self.type is None:
  if 1:
    if a_length<8: return False, "apical loop too short"
    if a_length>30: return False, "apical loop too long " #checking apical loop length
#  elif self.type ==1:
#    if a_length<8: return False, "apical loop (type 1) too short"
#    if a_length>23: return False, "apical loop (type 1) too long " #checking apical loop length for type 1
#  elif self.type ==2: 
#    if a_length<2: return False, "apical loop (type 2) too short"
#    if a_length>10: return False, "apical loop (type 2) too long " #checking apical loop length for type 2 
  
    
  s2_length=self.stem2_length()
  if s2_length<7:  return False, "stem2 too short"
  if s2_length>15: return False,  "stem2 too long" #checking stem2 length
  ins_s2=self.insertion_stem2()
  if max(ins_s2)-min(ins_s2) > 2:  return False, "too much bending" #too many more insertion on one side than the other
  if self.energy > opt['secis_energy']: return False, "free energy too high" #checking energy
  return True, "ok"

def default_filter(obj, min_conserved_per_side=3):
  """ obj can be a reverseblast or a SeB_exonerate"""
  if obj.distance_from_secis()>opt['max_secis_distance']:      return False
  if obj.category=='ALI_X' and obj.query_full_seq.count('X')>1: return False   
  if obj.protein().count('U')>1: return True
  pos_u_in_ali=obj.alignment.seq_of('t').index('U') #0-based
  conserved_left=0; conserved_right=0
  for i in range(pos_u_in_ali): 
    aa_query= obj.alignment.seq_of('q')[i];     aa_target= obj.alignment.seq_of('t')[i]
    if not '-' in  aa_query+aa_target and ( aa_query==aa_target or similar_aas(aa_query, aa_target) ): conserved_left+=1
  for i in range(pos_u_in_ali+1, obj.alignment.length()):
    aa_query= obj.alignment.seq_of('q')[i];     aa_target= obj.alignment.seq_of('t')[i]
    if not '-' in  aa_query+aa_target and (  aa_query==aa_target or similar_aas(aa_query, aa_target) ): conserved_right+=1
  if conserved_left< min_conserved_per_side: return False
  if conserved_right >= min_conserved_per_side: return True #####
  # conservation on the right is not wide. let's see if this fits the TRlike criteria
  if obj.query_prot_length - obj.query.boundaries()[1] <2 : #ultimate or penultimate residue of query is aligned
    try: 
      codon_downstream= upper(obj.downstream(0, 3).fasta_sequence()[1])
      if codon_downstream  in ['TAG', 'TAA']: return True  # and the codon just downstream is a non-TGA stop
    except: pass
  return False

def genes_in_the_middle_filter(b, candidate_blasth):
  """ to evaluate if a blast hit (maybe composed by multiple exons, superblasthit) can be accepted as a gene in the middle. it's important to keep this strict since it will filter out the original candidate. All the part where this is tested to be downstream of our candidate is done above, this is another test"""
  return b.evalue < 1e-10

def get_sequence_from_protein_db(title):
  global opt, prot_database
  write_to_file(title, temp_folder+'temp_gi')
  if not opt['net']:       
    bbash('fastacmd -d '+prot_database+' -i '+temp_folder+'temp_gi'+' > '+temp_folder+'temp_gi.fasta')
    full_title, seq = getfastalite(open(temp_folder+'temp_gi.fasta', 'r'))[0]
  else:              ### NEEDS TESTING! doens't work right now     
    fetched_hash = fetch_ncbi_sequences.main({'silent':1, 'f':1, 'm':'P', 'i':temp_folder+'temp_gi'})
    seq= fetched_hash[ fetched_hash.keys()[0] ]
  return full_title, seq 

def post_function_for_merging(upstream_g, downstream_g, out_g):
  """ used for blast hits"""
  out_g.query=out_g.query.boundaries_gene()
  out_g.evalue= min (upstream_g.evalue, downstream_g.evalue)
  out_g.cds_sequence=''
  if hasattr(upstream_g, 'is_candidate') or hasattr(downstream_g, 'is_candidate'): 
    out_g.is_candidate=True
    #checking if mergin removed a sec candidate
    if not 'u' in  out_g.protein() and ( 'u' in upstream_g.protein()+downstream_g.protein() ): printerr('WARNING candidate sec was removed by merging procedure in candidate number: '+out_g.id, 1)
  out_g.cds_sequence=None
  return out_g

def merge_by_colinearity_p2g_with_different_queries(p2g_list):
  """p2g_list is a  SeB_blast list """
  out_list=[]
  p2g_per_query={} #query: [p2gA, p2gB, ..]
  for p2g in p2g_list:
    p2g_per_query.setdefault( p2g.query.chromosome , [] )
    p2g_per_query[p2g.query.chromosome].append(p2g)

  for query_name in p2g_per_query:
    p2g_list_for_this_query = p2g_per_query[query_name]
    if len( p2g_list_for_this_query ) > 1:  merge_p2g_hits_by_colinearity(p2g_list_for_this_query, inplace=True, already_sorted=False, post_function=post_function_for_merging) #, get_sequence=get_sequence_from_protein_db) 
    out_list.extend(p2g_list_for_this_query)

  return out_list

####

class SeB_secis(infernalhit):
  def __init__(self, obj=None):
    """ Can be initialiased with a infernalhit or covelshit or secisearchhit"""
    gene.__init__(self)
    self.core=None; self.unpaired_nts=None; self.infernal_score=None; self.covels_score=None; self.secisearch_pattern=None; 
    self.alignment=alignment(); self.refined_ss=None; self.energy=None; self.energy_apical=None; self.energy_loop1=None;
    self.ss=None ; self.type= None
    self.cmali=infernal_model_alignment
    if not obj is None:
      if isinstance(obj, infernalhit):      self.load_infernal(obj)
      elif isinstance(obj, covelshit):      self.load_covels(obj)
      elif isinstance(obj, secisearchhit):  self.load_secisearch(obj)
      elif issubclass(obj.__class__, gene): self.load_gene(obj)
      else:  raise Exception, "SeB_secis-> load ERROR object not recognized: "+str(obj)

  def out_ss(self, aligned=False):
    """ returns a string with the secondary structure in output. It aligns to self.sequence (doesn't have the gaps in self.alignment), unless aligned=True"""
    o=''; gaps_target=0
    pairs=self.pairs()
    """
    for i, j, category in pairs:
      o+= '\n'+join( [ str(i), str(j), self.sequence()[i-1], self.sequence()[j-1], category], ' ')
    o+='\n target  '"""

    if self.refined_ss: ss= self.refined_ss
    else:               ss=self.ss
    bad_pairs=[p for p in pairs if p[2]=='o']
    ga_pairs=[p for p in pairs if p[2]=='k']
    left_of_bad_pairs={};     right_of_bad_pairs={}; left_of_ga_pairs={}; right_of_ga_pairs={}
    for p in bad_pairs:
      left_of_bad_pairs[p[0]]=True; right_of_bad_pairs[p[1]]=True
    for p in ga_pairs:
      left_of_ga_pairs[p[0]]=True;   right_of_ga_pairs[p[1]]=True
    for i in range(self.alignment.length()):
      if self.alignment.seq_of('t')[i]!='-':
        if i-gaps_target+1 in left_of_bad_pairs:    o+='{'
        elif i-gaps_target+1 in right_of_bad_pairs: o+='}'
        elif i-gaps_target+1 in left_of_ga_pairs:   o+='['
        elif i-gaps_target+1 in right_of_ga_pairs:  o+=']'
        else:                                     o+=ss[i]
      else:        
        gaps_target+=1
        if aligned: o+='-' 

    return o

  def load_infernal(self, hit):
    """ """
    self.chromosome= hit.chromosome; self.strand=hit.strand; self.id=hit.id;  self.target=hit.target
    self.exons=list(hit.exons) #hard copying exon list
    #filling extremes, then adding alignment
    hit.remove_Xs()
    seq_q=hit.alignment.seq_of('q'); seq_t=upper(hit.alignment.seq_of('t') )
    seq_q=nogap(self.cmali.rf)[   :hit.query.boundaries()[0]-1   ]+seq_q+  nogap(self.cmali.rf)[ hit.query.boundaries()[1]:  ]    
    seq_t=(hit.query.boundaries()[0]-1)*'-'+seq_t+len(nogap(self.cmali.rf)[ hit.query.boundaries()[1]:  ])*'-'
    ss_no_gaps=''
    for p_index, aa in enumerate(self.cmali.rf):
      if aa!='.': ss_no_gaps+=self.cmali.ss[p_index]
    ss=   ss_no_gaps[:hit.query.boundaries()[0]-1]+hit.ss+  ss_no_gaps[ hit.query.boundaries()[1]:  ]
    self.ss=''
    for c in ss:
      if c in '<[({':   self.ss+='('
      elif c in '>])}': self.ss+=')'
      else:             self.ss+='.'
    self.alignment.add('q', seq_q);     self.alignment.add('t', seq_t); 
    self.infernal_score=hit.score
    self.fill_core()
    if not opt['no_refine_ss']: 
      try: self.refine_structure()
      except: 
        write("WARNING can't refine a hit!", 1)
        self.refined_ss=None
    if opt['type']:        self.predict_stem3_with_RNAfold()
    self.compute_energy()
    self.fill_unpaired_nt()

  def load_covels(self, hit):
    """ Load covels hit, predict its secondary structure"""
    hit.extend(opt['covels_extension'], opt['covels_extension'], True)
    hit.check_boundaries(chromosome_lengths[hit.chromosome])
    self.load_gene(hit)
    self.covels_score=hit.score

  def load_gene(self, hit):
    """ Load from generic gene object, predict its secondary structure"""
    self.chromosome=hit.chromosome;  self.strand=hit.strand; self.target=hit.target
    self.exons=list(hit.exons)
    self.align_to_infernal_model(hit.fasta_sequence()[1])
    self.fill_core()
    if not opt['no_refine_ss']:
      try: self.refine_structure()
      except : 
        write("WARNING can't refine a hit!", 1)
        self.refined_ss=None
    if opt['type']:        self.predict_stem3_with_RNAfold()
    self.compute_energy()
    self.fill_unpaired_nt()

  def load_secisearch(self, hit):
    """ Load from secisearch object """
    reduction_left = len(hit.seq.split()[0] )
    reduction_right= len(hit.seq.split()[-1])
    reduced_seq    = join(hit.seq.split()[1:-1], '')
    self.chromosome=hit.chromosome;  self.strand=hit.strand; self.target=hit.target
    self.exons=[];     self.add_exon( hit[0][0]+reduction_left, hit[0][1]-reduction_right)
    self.secisearch_pattern=hit.pattern
    self.align_to_infernal_model(reduced_seq)
    self.fill_core()
    if not opt['no_refine_ss']:
      try: self.refine_structure()
      except: 
        write("WARNING can't refine a hit!", 1)
        self.refined_ss=None
    if opt['type']:        self.predict_stem3_with_RNAfold()      
    self.compute_energy()
    self.fill_unpaired_nt()

  def align_to_infernal_model(self, sequence):
    """ Uses the global variable infernal_cm, it aligns the sequence (self.seq if not provided as argument) to the infernal cm model using cmalign, and fills .alignment, and .ss using the information in cmalign out; used for covels and secisearch predictions """
    temp_fasta_for_cmalign=temp_folder+'cmalign_inputfile.fa';temp_out_cmalign=temp_folder+'cmalign_outfile'
    write_to_file(">hit\n"+sequence, temp_fasta_for_cmalign)
    bbash('cmalign '+infernal_cm+' '+temp_fasta_for_cmalign+' > '+temp_out_cmalign)
    a=alignment(); a.load_cmalign_out(temp_out_cmalign)
    self.alignment.add('q', replace(a.rf, '.', '-')); self.alignment.add('t', a.seq_of('hit'))
    self.ss=''
    for c in a.ss: 
      if c in '<[({':   self.ss+='('
      elif c in '>])}': self.ss+=')'
      else:             self.ss+='.'

  def fill_core(self):
    """ fill the .core and .unparied_nts with positions (see below for meaning) """
    self.core=[[], []]  ## the UGAz - zGAU positions , in form [ [posU, posG, posA, posZ], [posZ, posG, posA, posU]  ] ; if not all of them are aligned, anyway the first and second list are same length
    pairs=self.infernal_pairs()
    ### compute core   #all pos are 1-based!
    for i, j, category in pairs:
      i_in_ali=  self.alignment.position_in_ali('t', i)
      i_in_query=self.alignment.position_in_seq('q', i_in_ali)
      if i_in_query in [17, 18, 19, 20]: 
        self.core[0].append(i); self.core[1].insert(0, j)
    
    if not any( self.core ): self.core=None

  def fill_unpaired_nt(self):
    """  compute unpaired_nts, 1 based positions (2 or 3 nt in total) _ > NB using refined ss if available"""
    try: 
      apical_loop_bounds=self.apical_loop_boundaries()
      if not apical_loop_bounds: return None
      portion=upper(self.sequence()[ apical_loop_bounds[0]-1: apical_loop_bounds[1]  ]) [:5] #allowing max 2 nts before the unpaired nts
      self.unpaired_nts=None
      for nt in ('AAA', 'AA', 'CC'): #trying to find AAA, or AA, or CC. stopping at the first one found
        find_in_portion=portion.find(nt)
        if find_in_portion!=-1 and ( (find_in_portion<=3 and self.apical_loop_length()>=9) or (find_in_portion<=2) ): # limit: 3 before unpaired if loop is more than 8, 2 otherwise
          self.unpaired_nts=  [ p+find_in_portion+apical_loop_bounds[0]     for p in range(len(nt))  ]
          break        
    except: 
      printerr('WARNING fill_unpaired_nt failed on a SECIS prediction!', 1)

  def apical_loop_length(self, extend_to_mismatches=True):
    """ using apical_loop_boundaries"""
    al_boundaries=self.apical_loop_boundaries(extend_to_mismatches=True)
    if not al_boundaries: return 0
    return al_boundaries[1]-al_boundaries[0]+1

  def apical_loop_boundaries(self, extend_to_mismatches=True):
    """ return [start, stop] or returns None if no apical loop can be identified (there's no stem2). If the structure has been refined (by default it happens), extend to mismatches makes no sense """
    s2_pairs=self.stem2_pairs()
    if s2_pairs:
      i=1 #to count as apical loop also the bad pairs at the end of stem2
      if extend_to_mismatches:
        while s2_pairs[-i][2]=='o': i+=1
      a= s2_pairs[-i][0]
      b= s2_pairs[-i][1]
      if b-a > 1: return a+1, b-1
    return None
  
  def loop1_boundaries(self, extend_to_mismatches=True):
    """ Computed on the secondary structure by infernal, to facilitate the calling. returns None if nothing is aligned to that portion. otherwise, a list of [ [pos5_start, pos5_end], [pos3_start, pos3_end]] is returned, all 1 based"""
    s2_pairs=self.stem2_pairs()
    if not s2_pairs or not self.core: return None
    ip=self.infernal_pairs()
    possible_extension_5=[]; possible_extension_3=[]
    last_stem1_pairs_i=None
    for i, j, category in ip:
      i_in_query= self.alignment.position_in_seq('q',    self.alignment.position_in_ali('t', i))
      if i_in_query>=17: break
      if category!='o' or not extend_to_mismatches: 
        possible_extension_5=[]; possible_extension_3=[]
      else: 
        possible_extension_5.append(i)
        possible_extension_3.insert(0,j)
      last_stem1_pairs_i=i;       last_stem1_pairs_j=j
    if last_stem1_pairs_i:
      #now I have in last_stem1_pairs_i the position of the last pair of stem1, and in possible_extension the consecutive positions with bad pairing starting from this, backwards.
      loop_positions_5=possible_extension_5
      loop_positions_5.extend(   range(last_stem1_pairs_i+1, min(self.core[0]))   )   #extending with positions of loop1, 5' portion
      loop_positions_3=[]
      loop_positions_3=range( max(self.core[1])+1, last_stem1_pairs_j )
      loop_positions_3.extend(possible_extension_3)
      return [ [min(loop_positions_5), max(loop_positions_5)], [min(loop_positions_3), max(loop_positions_3)] ]

  def stem2_length(self):
    """Computed on the secondary structure by infernal, to facilitate the calling of what's stem2. Not counting the core as part of stem2 here! """
    pairs_on_stem2= self.stem2_pairs()
    return len(  [ 1  for obj in pairs_on_stem2 if obj[2]!='o' ]  )

  def stem2_pairs(self):
    out=[]
    for i, j, category in self.pairs(from_infernal=self.refined_ss is None):
      i_in_query= self.alignment.position_in_seq('q',    self.alignment.position_in_ali('t', i))
      if i_in_query>19 and i_in_query<33 and (self.core is None or not i in self.core[0]):        out.append([i, j, category])
    return out

  def stem1_pairs(self):
    out=[]
    for i, j, category in self.pairs(from_infernal=self.refined_ss is None):
      i_in_query= self.alignment.position_in_seq('q',    self.alignment.position_in_ali('t', i))
      if i_in_query<19 and (self.core is None or not i in self.core[0]):        out.append([i, j, category])
    return out

  def insertion_stem2(self):
    """ Returns [x, y] where x is the number of insertion in the first part (5') of the stem 2, and y is the number of insertion in the 3' part.  """
    pairs_on_stem2= self.stem2_pairs()
    if not pairs_on_stem2: return [0 , 0]
    insertions=self.insertions()
    insertions_left=[]
    insertions_right=[]
    for i in insertions:
      if i >= pairs_on_stem2[0][0]-1 and i< pairs_on_stem2[-1][0]: insertions_left.append(i) #including the position just before stem2, as we should lcount them as insertions
      elif i>pairs_on_stem2[-1][1] and i<=pairs_on_stem2[0][1]+1:  insertions_right.append(i) #idem  
    return [len(insertions_left), len(insertions_right)]


  def fasta_with_corrected_infernal_ss(self):
    """ return a fasta like string with sequence in target and secondary structure. The secondary structure is corrected to have all matches -> equal number of (  and ) """
    a=alignment()
    ss=self.ss
    left_opened_parenthesis_gap=0;     right_opened_parenthesis_gap=0;     left_closed_parenthesis_gap=0;     right_closed_parenthesis_gap=0
    pre_seq=upper(self.alignment.seq_of('t'))
    tot_pairs_in_model=ss.count('(')
    index_open_par=0; index_close_par=tot_pairs_in_model+1; must_unpair={}
    for index, char in enumerate(ss):
      if ss[index]=='(':
        index_open_par+=1
        if pre_seq[index]=='-':
          must_unpair[index_open_par]=True
      if ss[index]==')':
        index_close_par-=1
        if pre_seq[index]=='-':
          must_unpair[index_close_par]=True
    pre_ss=ss #we will not change pre_ss, while ss will be modified to match number of ( and )
    already_unpaired=0
    for index_pair in sorted(must_unpair.keys()): # 1 means the first pair (first on 5' with last on 3')
      pos_open_to_unpair=  pre_ss.replace('(', '&', index_pair-1).find('(') #pos is 0 based
      pos_close_to_unpair= len(pre_ss) -1 - pre_ss[::-1].replace(')', '&', index_pair-1).find(')')  #reversing ss to make the replace-find trick work
      ss=ss[:pos_open_to_unpair]+'.'+ss[pos_open_to_unpair+1:pos_close_to_unpair]+'.'+ss[pos_close_to_unpair+1:]
    out_seq=''; out_ss=''
    for index, char in enumerate(pre_seq): 
      if char!='-':
        out_seq+=char
        out_ss+=ss[index]
    return ">"+self.header()+'\n'+out_seq+'\n'+out_ss

  def RNAplot_options(self):
    """ Return the adds to the commandline to put for getting a nice SECIS image. must be added to: -t 1  < inputfile  """

    #some of this code if left from the phase in which I was testing colors
                            #cyan                      #cyan more blue
    color_per_category={'a':'0.492', 'c':'0.492', 'w':'0.492', 'k':'0.583', 'o':'0.00'} 
                       #green
    color_sat_core=' 0.333 1 ' # if color_sat_core is not False, the core will be displayed with this color. if it's false, it will be displayed with the same colors as color_per_category dictates, but with the saturation in saturation_core
    color_sat_nts='  0.333 1 ' ;     saturation_core=' 1 ';     saturation_rest=' 0.28 '
    o=' --pre "'
    pairs=self.pairs()
    ## coloring core
    if color_sat_core and  self.core:
      for index_core in range(len(self.core[0])):
        i=self.core[0][index_core]; j=self.core[1][-1-index_core]
        o+=' '+str(i)+' '+str(j)+color_sat_core+' colorpair '
      if len(self.core[0])==4: 
        if self.core[0][0]!=1:
          o+=' '+str(self.core[0][0]-1)+' '+str(self.core[0][0]-1)+color_sat_core+' colorpair ' #base preceding the paired core
    for i, j, category in pairs:
      if  self.core and (i in self.core[0] or j in self.core[1] or i==self.core[0][0]-1 ): #part of the core or preceding base
        if not color_sat_core:        o+=' '+str(i)+' '+str(j)+' '+color_per_category[category]+' '+saturation_core+' colorpair '  #saturation for core ## drawing only if color schems wants same color for core but more saturation ###NOT USED RIGHT NOW
      else:                  
        if category=='o':
          #coloring mismatch as two separate bubbles
          o+=' '+str(i)+' '+str(i)+' '+color_per_category[category]+' '+saturation_rest+' colorpair '  #mismatch
          o+=' '+str(j)+' '+str(j)+' '+color_per_category[category]+' '+saturation_rest+' colorpair '  #mismatch
        else:
          o+=' '+str(i)+' '+str(j)+' '+color_per_category[category]+' '+saturation_rest+' colorpair '  #saturation for other pairs
    if self.unpaired_nts:
      for p in self.unpaired_nts:
        o+=' '+str(p)+' '+str(p)+color_sat_nts+' colorpair '
        o+=' '+str(p)+' cmark ' # circles around unpaired nts
    for ins_pos in self.insertions():
      o+=' '+str(ins_pos)+' '+str(ins_pos)+' '+color_per_category['o']+' '+saturation_rest+' colorpair '  #insertion as mismatch
        
    o+='"'
    return o

  def infernal_pairs(self):
    return self.pairs(from_infernal=True)

  def pairs(self, from_infernal=False):
    """ returns a list of the nucleotides paired in the infernal structure, and their "category"   ; output is a list of [index1, index2, "X"], where X can be: 
    a -> AT pair    c -> CG pair   k -> GA pair w-> GU pairs   o-> others
    indexes in output are 1-based
    """
    if from_infernal or self.refined_ss is None:       a,  seq, ss = self.fasta_with_corrected_infernal_ss().split('\n')
    else:                                              seq=self.sequence(); ss=replace(self.refined_ss, '-', '')

    out=[]    
    closed_par_already_matched={} #1-based
    last_found=   ss.rfind('(')
    while last_found!=-1:
      i_search=0
      while ss[last_found+i_search]!=')' or closed_par_already_matched.has_key(last_found+i_search+1): i_search+=1        #scanning for pairing parenthesis
      x=last_found+1;     y=last_found+i_search+1 #pos 1 based of matchin pairs: x-y is a pair, with x<y
      closed_par_already_matched[y]=1
      try: category=category_per_pair.get(seq[x-1], seq[y-1] )
      except KeyError:  category='o'
      out.append([x, y, category])
      last_found=   ss.rfind('(', 0, last_found)
    out.sort(key=lambda x:x[0])

    return out


  def refine_structure(self):
    """ Refine the structure by infernal by correcting most obvious mistakes"""
    t, seq, ss= self.fasta_with_corrected_infernal_ss().split('\n')
    apical_bounds=self.apical_loop_boundaries(False)
    loop1_bounds=self.loop1_boundaries(False)
    pairs=self.infernal_pairs()

    if loop1_bounds:
      #counting bad pairs at the end of stem1
      bad_pairs_at_end_of_stem1=[]
      for i, j, category in pairs:
        if i>=loop1_bounds[0][0]: break
        if category in 'ok': bad_pairs_at_end_of_stem1.append([i,j])
        else: bad_pairs_at_end_of_stem1=[]
      #extending loop1 up
      extension_loop1_up=[] #positions list [i, j] , 1 based
      if  not len(bad_pairs_at_end_of_stem1)>=2:
        i_index=loop1_bounds[0][0]
        j_index=loop1_bounds[1][1]
        while i_index <=loop1_bounds[0][1] and (self.core is None or (not i_index in self.core[0] and not j_index in self.core[1] )):
          i_nt=self.sequence()[i_index-1]
          j_nt=self.sequence()[j_index-1]
          try: category=category_per_pair.get(i_nt, j_nt)
          except KeyError:  category='o'      
          if category in 'acw': 
            extension_loop1_up.append([i_index, j_index])
            i_index+=1;           j_index-=1
          else:          break
        for i, j in extension_loop1_up:        ss=ss[:i-1]+'('+ss[i:j-1]+')'+ss[j:]
      if not extension_loop1_up and  bad_pairs_at_end_of_stem1: #if no extension was done and there are bad pairs at the end, unpairing them.
        for i, j in bad_pairs_at_end_of_stem1:
          ss=ss[:i-1]+'.'+ss[i:j-1]+'.'+ss[j:]

      #coutning bad paris at start of stem1
      bad_pairs_at_start_of_stem1=[]
      for i, j, category in pairs:
        if i>=loop1_bounds[0][0]: break
        if category in 'ok': bad_pairs_at_start_of_stem1.append([i,j])
        else: break    
      #extending loop1 down  
      extension_loop1_down=[]
      if not len(bad_pairs_at_start_of_stem1)>=2:
        i_index=pairs[0][0]-1
        j_index=pairs[0][1]+1
        while i_index>=1 and j_index<=len(self.sequence()):
          i_nt=self.sequence()[i_index-1]
          j_nt=self.sequence()[j_index-1]
          try: category=category_per_pair.get(i_nt, j_nt)
          except KeyError:  category='o'
          if category in 'acw':
              extension_loop1_down.append([i_index, j_index])
              i_index-=1;           j_index+=1
          else:          break
        for i, j in extension_loop1_down:        ss=ss[:i-1]+'('+ss[i:j-1]+')'+ss[j:]    
      if not extension_loop1_down and  bad_pairs_at_start_of_stem1:
        for i, j in bad_pairs_at_start_of_stem1:
          ss=ss[:i-1]+'.'+ss[i:j-1]+'.'+ss[j:]

    if apical_bounds:
      #extending apical loop
      bad_pairs_at_end_of_stem2=[]
      for i, j, category in pairs:
        if i>=apical_bounds[0]: break
        if category in 'ok': bad_pairs_at_end_of_stem2.append([i,j])
        else: bad_pairs_at_end_of_stem2=[]
      extension_apical_loop=[] #positions list [i, j] , 1 based  
      if not len(bad_pairs_at_end_of_stem2)>0: ## differnet number! 1 as threshold
        i_index=apical_bounds[0]
        j_index=apical_bounds[1]
        while i_index <=apical_bounds[1] and i_index<j_index:
          i_nt=self.sequence()[i_index-1]
          j_nt=self.sequence()[j_index-1]
          try: category=category_per_pair.get(i_nt, j_nt)
          except KeyError:  category='o'
          if category in 'acw':
            extension_apical_loop.append([i_index, j_index])
            i_index+=1;           j_index-=1
          else:          break
        for i, j in extension_apical_loop:        ss=ss[:i-1]+'('+ss[i:j-1]+')'+ss[j:]
      if not extension_apical_loop and  bad_pairs_at_end_of_stem2: #if no extension was done and there are bad pairs at the end, unpairing them.
        for i, j in bad_pairs_at_end_of_stem2:
          ss=ss[:i-1]+'.'+ss[i:j-1]+'.'+ss[j:]

    a=alignment(); a.add('SS', ss); a.add('t', seq);     b=a.transfer_alignment( self.alignment , dont_shrink=True )
    self.refined_ss=b.seq_of('SS') 
    
    #refining insertions
    self.refine_insertions()

  def refine_insertions(self, starting_point=1):
    """ call only after refine_structure. It fixes insertions by trying to match with another inserted nt on the other side of the stem, or by unpairing a bad pair which is close to the insertion either above or below the insertion.   Basically it makes more logical the predictions by infernal"""
    self.old_refined_ss=self.refined_ss 
    ss=replace(self.refined_ss, '-', '')
    def repeat_function(starting_point=1):
      close()
      self.refine_insertions(starting_point=starting_point)
    def close():
      a=alignment(); a.add('SS', ss); a.add('t', self.sequence());     b=a.transfer_alignment( self.alignment , dont_shrink=True )
      self.refined_ss=b.seq_of('SS')

    refined_pairs=self.pairs()
    insertions=self.insertions()
    index_pair=0 #0-based
    ins_index=0
    while ins_index < len(insertions):
      ins=insertions[ins_index]
      if ins<starting_point: 
        ins_index+=1
        continue
      ins_is_left= ins <= refined_pairs[-1][0]
      if ins_is_left:
        while refined_pairs[index_pair][0]+1<ins: index_pair+=1
        #now refined_pairs[index_pair] is the pair just before the insertion
        length_ins_right=0
        while refined_pairs[index_pair][1]-length_ins_right-1 in insertions: length_ins_right+=1
        length_ins_left=1
        while ins+length_ins_left in insertions: length_ins_left+=1

        if length_ins_right:
          #we have an insertion in the same position on the right side of the stem.
          best_match_pair=None
          for i in range(ins, ins+length_ins_left):
            for j_insertion_index in range(length_ins_right):
              j= refined_pairs[index_pair][1]-j_insertion_index-1
              i_nt=self.sequence()[i-1];     j_nt=self.sequence()[j-1]
              try: category=category_per_pair.get(i_nt, j_nt)
              except KeyError:  category='o'
              if best_match_pair is None or category_per_pair_hierarchy.index(category) < category_per_pair_hierarchy.index(best_match_pair[2]):
                best_match_pair=[i, j, category]

          #matching the best pair in the bulge. we don't expect more than one
          i, j, category=best_match_pair
          ss=ss[:i-1]+'('+ss[i:j-1]+')'+ss[j:]
          repeat_function(ins+1)
          return
        #checking the two close pairs to see if there are bad pairs that can be fixed
        else:

          category_gain_below=0
          if refined_pairs[index_pair][0]==ins-1:
          #now trying to pair this left insertion with pair on the right just below
            i_below=ins;      j_below=refined_pairs[index_pair][1]
            i_nt_below=self.sequence()[i_below-1];     j_nt_below=self.sequence()[j_below-1]
            try: category_below=category_per_pair.get(i_nt_below, j_nt_below)
            except KeyError:  category_below='o'
            category_gain_below=category_per_pair_hierarchy.index( refined_pairs[index_pair][2]  ) - category_per_pair_hierarchy.index(category_below) 

          category_gain_above=0
          if index_pair+1<len(refined_pairs) and refined_pairs[index_pair+1][0]==ins+1: ## now trying with the pair just above
            i_above=ins;      j_above=refined_pairs[index_pair+1][1]
            i_nt_above=self.sequence()[i_above-1];     j_nt_above=self.sequence()[j_above-1]
            try: category_above=category_per_pair.get(i_nt_above, j_nt_above)
            except KeyError:  category_above='o'
            category_gain_above=category_per_pair_hierarchy.index( refined_pairs[index_pair][2]  ) - category_per_pair_hierarchy.index(category_above) 
          
          if category_gain_below>0 and category_gain_below>=category_gain_above:
            #improving using pair below
            ss=ss[:refined_pairs[index_pair][0]  -1]+'.'+ss[refined_pairs[index_pair][0]:refined_pairs[index_pair][1]-1]+'.'+ss[refined_pairs[index_pair][1]:]     #unpair the old one
            ss=ss[:i_below-1]+'('+ss[i_below:j_below-1]+')'+ss[j_below:] #pair the new one
            repeat_function(ins+1)
            return
          elif category_gain_above >0:
            #improving using pair above
            ss=ss[:refined_pairs[index_pair+1][0]  -1]+'.'+ss[refined_pairs[index_pair+1][0]:refined_pairs[index_pair+1][1]-1]+'.'+ss[refined_pairs[index_pair+1][1]:]     #unpair the old one
            ss=ss[:i_above-1]+'('+ss[i_above:j_above-1]+')'+ss[j_above:] #pair the new one
            repeat_function(ins+1)          
            return
          else: pass #not improving. leaving things as they are  

      else:
          #insertion is on the right. bulges are wlready fixed, but we have to check if we can improve bad pairs
          index_pair=0
          while refined_pairs[index_pair][1]-1>ins: index_pair+=1
          #now refined_pairs[index_pair] is the pair just before (below) the insertion          


          category_gain_below=0
          if refined_pairs[index_pair][1]==ins+1:
          #now trying to pair this left insertion with pair on the right just below
            i_below=refined_pairs[index_pair][0];    j_below=ins
            i_nt_below=self.sequence()[i_below-1];     j_nt_below=self.sequence()[j_below-1]
            try: category_below=category_per_pair.get(i_nt_below, j_nt_below)
            except KeyError:  category_below='o'
            category_gain_below=category_per_pair_hierarchy.index( refined_pairs[index_pair][2]  ) - category_per_pair_hierarchy.index(category_below) 

          category_gain_above=0
          if index_pair+1<len(refined_pairs) and refined_pairs[index_pair+1][1]==ins-1: ## now trying with the pair just above
            i_above=refined_pairs[index_pair+1][0];    j_above=ins 
            i_nt_above=self.sequence()[i_above-1];     j_nt_above=self.sequence()[j_above-1]
            try: category_above=category_per_pair.get(i_nt_above, j_nt_above)
            except KeyError:  category_above='o'
            category_gain_above=category_per_pair_hierarchy.index( refined_pairs[index_pair][2]  ) - category_per_pair_hierarchy.index(category_above) 
          

          if category_gain_below>0 and category_gain_below>=category_gain_above:
            #improving using pair below
            ss=ss[:refined_pairs[index_pair][0]  -1]+'.'+ss[refined_pairs[index_pair][0]:refined_pairs[index_pair][1]-1]+'.'+ss[refined_pairs[index_pair][1]:]     #unpair the old one
            ss=ss[:i_below-1]+'('+ss[i_below:j_below-1]+')'+ss[j_below:] #pair the new one
            repeat_function(ins+1)
            return
          elif category_gain_above >0:
            #improving using pair above
            ss=ss[:refined_pairs[index_pair+1][0]  -1]+'.'+ss[refined_pairs[index_pair+1][0]:refined_pairs[index_pair+1][1]-1]+'.'+ss[refined_pairs[index_pair+1][1]:]     #unpair the old one
            ss=ss[:i_above-1]+'('+ss[i_above:j_above-1]+')'+ss[j_above:] #pair the new one
            repeat_function(ins+1)          
            return
          else: pass #not improving. leaving things as they are  

      ins_index+=1
    close()

  def insertions(self):
    """returns a list of positions 1 based of insertions """
    out=[]
    s1_pairs=self.stem1_pairs()
    if s1_pairs:
      s1_pairs_left= [p[0] for p in s1_pairs]
      s1_pairs_right=[p[1] for p in s1_pairs]
      insertions_left_s1=  [ i   for i in range(s1_pairs[0][0], s1_pairs[-1][0]+1)     if not i in s1_pairs_left]
      insertions_right_s1= [ j   for j in range(s1_pairs[-1][1], s1_pairs[0][1]+1)     if not j in s1_pairs_right]  
      out.extend(insertions_left_s1)
      out.extend(insertions_right_s1)
    s2_pairs=self.stem2_pairs()
    if s2_pairs:      
      s2_pairs_left= [p[0] for p in s2_pairs]
      s2_pairs_right=[p[1] for p in s2_pairs]
      insertions_left_s2=  [ i   for i in range(s2_pairs[0][0], s2_pairs[-1][0]+1)     if not i in s2_pairs_left]
      insertions_right_s2= [ j   for j in range(s2_pairs[-1][1], s2_pairs[0][1]+1)     if not j in s2_pairs_right]
      out.extend(insertions_left_s2)
      out.extend(insertions_right_s2)
    return out
    
    
  def compute_energy(self):
    """ with RNAeval """
    if self.refined_ss:      ss_no_gaps= replace(self.refined_ss, '-', '')
    else:                    
      t, seq, ss_no_gaps= self.fasta_with_corrected_infernal_ss().split('\n')
    #computing free energy
    out_rnaeval=bbash('echo "'+self.sequence()+'\n'+ss_no_gaps+'" | RNAeval ')
    self.energy=float(out_rnaeval.split('\n')[-1].split()[-1].strip('()'))

  def predict_stem3_with_RNAfold(self):
    """ #keeping this just in case we get back to it
    t, seq, ss= self.fasta_with_corrected_infernal_ss().split('\n')
    pairs=self.infernal_pairs()
    #mismatching non canonical-or-wooble pairs ; 
    for i, j, category in pairs:
      if category=='o' and (self.core is None or not i in self.core[0]):        ss=ss[:i-1]+'.'+ss[i:j-1]+'.'+ss[j:]
    #forcing apicel  unpaired nucleotides to be unpaired
    if self.unpaired_nts:
      for p in self.unpaired_nts:        ss=ss[:p-1]+'x'+ss[p:] 
    """
    
    #refining secondary strucure of apical loop
    if self.refined_ss is None:       a,  seq, ss = self.fasta_with_corrected_infernal_ss().split('\n') # only ss is used
    else:                             ss=replace(self.refined_ss, '-', '')

    a= self.apical_loop_boundaries()
    refined_ss_no_gaps=ss
    if a: 
      apical_start, apical_end= a[0], a[1]
      seq_apical=self.sequence()[apical_start-1:apical_end]      
      out_rnafold=bbash('echo "'+'xxx'+seq_apical+'xxx" | RNAfold --noPS ')
      ss_apical_loop=out_rnafold.split('\n')[-1].split()[0][3:-3]

      if ss_apical_loop.count('(')<2: 
        self.type=1
        return    ### NB: keeping the structure only if at least two nucleotides are paired
      
      self.type=2
      self.energy_apical= float(out_rnafold.split('\n')[-1].split()[-1].strip('()'))
      refined_ss_no_gaps=refined_ss_no_gaps[:apical_start-1]+ss_apical_loop+refined_ss_no_gaps[apical_end:]


      a=alignment(); a.add('SS', refined_ss_no_gaps); a.add('t', self.sequence());     b=a.transfer_alignment( self.alignment , dont_shrink=True )
      self.refined_ss=b.seq_of('SS')

    """ ### old code, not using this: see refine_structure instead.
    #refining secondary structure of loop1  
    loop1_boundaries=self.loop1_boundaries()
    if loop1_boundaries:
      seq_loop_5=self.sequence()[loop1_boundaries[0][0]-1:loop1_boundaries[0][1]]
      seq_loop_3=self.sequence()[loop1_boundaries[1][0]-1:loop1_boundaries[1][1]]
      out_rnafold=bbash('echo "'+'xxx'+seq_loop_5+'xxxxx'+seq_loop_3+'xxx" | RNAfold --noPS ')
      self.energy_loop1= float(out_rnafold.split('\n')[-1].split()[-1].strip('()'))
      ss_loop1=out_rnafold.split('\n')[-1].split()[0][3:-3] #careful! contains central xxxxx still
      refined_ss_no_gaps= refined_ss_no_gaps[:loop1_boundaries[0][0]-1]+ ss_loop1[:len(seq_loop_5)]+refined_ss_no_gaps[loop1_boundaries[0][1]:loop1_boundaries[1][0]-1]+ss_loop1[len(seq_loop_5)+5:]+refined_ss_no_gaps[loop1_boundaries[1][1]:]
    """ 

  def draw_image(self, to_file):
    """ Draw a png image into a file """
    for program in ['RNAplot', 'convert']:
      if 'not found' in bash(program+' -h')[1]: raise Exception, "SeB_secis->draw_image ERROR can't draw image: program "+program+" not found "
    ## if a refined structure is present, drawing that, otherwise, drawing infernal ss
    if not self.refined_ss is None:          write_to_file(">"+self.id+'\n'+self.sequence()+'\n'+replace(self.refined_ss, '-', '') , temp_folder+'drawing_secis.fad'  )
    else:                                write_to_file(self.fasta_with_corrected_infernal_ss(), temp_folder+'drawing_secis.fad')
    bbash('> '+to_file) #flushing output file
    ps_filename=temp_folder+str(self.id)+'_ss.ps'
    layout_flag='1' if not opt['fd'] else '0'
    bbash('cd '+temp_folder+'; RNAplot -t '+layout_flag+' < '+'drawing_secis.fad '+self.RNAplot_options()) #this will create a file called $ps_filename ### -> colors in RNAplot_options()
    transformed_ps_filename=temp_folder+str(self.id)+'_ss.transformed.ps'
    write_to_file(  transform_secis_ps(ps_filename), transformed_ps_filename )
    bbash('convert -rotate 180  -density '+str(opt['dpi_resolution'])+' -geometry 100%  '+transformed_ps_filename+' '+to_file)   

  def grade(self):
    """Assign a grade (categorized score) which goes from C to A, where A is the best. This depends on the presence of unpaired AA, bending coefficient, mismatches in stem2 """
    score=0.0
    i=self.insertion_stem2()
    pairs_stem2=self.stem2_pairs()
    if  max(i)-min(i)<=1: score+=1           #bending coefficient 1 or 0 for stem2
    if  len([1 for obj in pairs_stem2 if obj[2]=='o']) < 3:      score+=1 #bad pairs on stem2 <3 
    if self.unpaired_nts:
      nts=join([ self.sequence()[p-1]  for p in self.unpaired_nts ], '')
      if 'AA' in nts:             score+=1       # AA or AAA in unpaired nts
      elif 'CC' in nts:           score+=0.4     # CC in unpaired nts
    if self.get_cove_score() >= 15:     score+=1 # great cove score
    ## now penalizing some features
    if self.core: 
      #checking position after core
      last_core_pos_x=max(self.core[0]);    last_core_pos_y=min(self.core[1])
      after_pos_x_is_paired=False; after_pos_y_is_paired=False;
      for x, y, category in pairs_stem2:
        if x == last_core_pos_x+1:
          after_pos_x_is_paired=True
          if y == last_core_pos_y-1: after_pos_y_is_paired=True
          if category=='o':           score-=0.6 # mismatch after core
          break
      if not after_pos_x_is_paired or not after_pos_y_is_paired:     score-=0.6 #insertion after core on either side
       #checking position before core
      first_core_pos_x=min(self.core[0]);   first_core_pos_y=max(self.core[1]);
      if first_core_pos_x>1 and len(self.sequence())>first_core_pos_y:
        nt_precore_x= self.sequence()[first_core_pos_x-2]
        nt_precore_y= self.sequence()[first_core_pos_y]
        try: category=category_per_pair.get(nt_precore_x, nt_precore_y)
        except KeyError:  category='o'
        if category in 'acw':    score-=0.4 #match just before core

    last_pair_category=''
    for x, y, category in pairs_stem2:
      if last_pair_category and last_pair_category=='o' and category=='o':   score-=1  #two bad pairs are close (O-filter) 
      last_pair_category=category

    if score<1.5:    out='C'
    elif score<2.5:  out='B'
    elif score>=2.5 and self.get_cove_score() >= 15:         out='A'
    else: out='B'
    return out

  def get_cove_score(self):
    """ compute a covels score is not computed already"""
    if self.covels_score is None:
      tempfileseq=  temp_folder+'temp_to_compute_cove.fa'
      tempcovelsout=temp_folder+'temp_to_compute_cove.coveout'
      write_to_file('>'+self.id+'\n'+self.sequence(), tempfileseq)
      bbash('covels-SE '+covels_cm+' '+tempfileseq+' > '+tempcovelsout)
      covels_hits=sorted( parse_covels(tempcovelsout).all(), key=lambda x:x.score, reverse=True  ) #in almost all cases we'll have a single covels hit, or none. anyway i script the general case
      if covels_hits:
        best_covels_hit = covels_hits[0]
        self.covels_score = best_covels_hit.score
      else: self.covels_score=0.0
    return self.covels_score

  def summary(self, add=None, with_grade=True):
    out= '------- SECIS element -- id:'+str(self.id)
    if with_grade: out+=' -- grade:'+str(self.grade())
    out+=' ----\n'
    out+=' Chromosome = '+self.chromosome
    out+=' Strand = '+self.strand+' Positions = '
    if self.strand=='+': out+=str(self.boundaries()[0])+' - '+str(self.boundaries()[1])+'\n'
    else:                out+=str(self.boundaries()[1])+' - '+str(self.boundaries()[0])+'\n'
    out+=' Found by: '
    if not self.infernal_score is None:     out+=" Infernal (score="+str(self.infernal_score)+")"
    if not self.covels_score is None and self.covels_score>0:       out+=" Covels (score="+str(self.covels_score)+")"
    if not self.secisearch_pattern is None: out+=" SECISearch1.0 (pattern="+str(self.secisearch_pattern)+")"
    out+="\n"
    out+=" Stem2 length = "+str(self.stem2_length())+"   Apical loop length = "+str(self.apical_loop_length())+'   Insertion in stem2 = '+str(self.insertion_stem2())[1:-1]+'\n'
    if not self.energy is None:             
      out+=" Free Energy = "+str(self.energy)
      try: # this is to deal with SeB objects computed before adding .energy* attribute
        if not self.energy_apical is None:        out+="   Free energy of apical portion = "+str(self.energy_apical)
        if not self.energy_loop1 is None:        out+="   Free energy of loop1  = "+str(self.energy_loop1)
      except: pass
      out+='\n'
    try:   # this is to deal with SeB objects computed before adding .type attribute
      if not self.type is None:    out+=" Type = "+str(self.type)+'\n'
    except: pass
      
    if add: out+=add+'\n'
    out+='\n'
    out+=self.pretty_alignment()
    return out

  def pretty_alignment(self):
    out=''
    if self.refined_ss:     out+='         '+self.ss+'\n'  #if two ss predictions are available, we show for one of the model above the model
    out+=' model   '+self.alignment.seq_of('q')+'\n'
    out+='         '+self.out_ss(aligned=True)+'\n'
    out+=' target  '+self.alignment.seq_of('t')+'\n'
    out+='         '+self.annotation_line()
    return out

  def annotation_line(self):
    """ Returns a string of same length of the seq""" 
    a_line=''
    positions_core=[]; positions_unpaired_nt=[]
    if self.core:      positions_core=self.core[0]+self.core[1]
    if self.unpaired_nts: positions_unpaired_nt=self.unpaired_nts
    gaps_target=0
    for pos in range(self.alignment.length()):
      if self.alignment.seq_of('t')[pos]=='-': 
        gaps_target+=1
        a_line+=' '
      else:
        if   pos-gaps_target+1 in positions_core:          a_line+='*'  #core position
        elif pos-gaps_target+1 in positions_unpaired_nt:   a_line+='^'  #unpaired nt position     
        else:                                              a_line+=' '                             
    return a_line

  def fasta(self):
    """ Returns the text for a fasta file"""
    header=">SECIS."+self.header(no_target=True)
    if not self.type is None: header+=" Type="+str(self.type)
    return header+'\n'+self.sequence()

  def gff(self):
    """ Returns the text for a gff file"""    
    comment=''
    if not self.infernal_score is None:     comment+="Infernal_score="+str(self.infernal_score)+" "
    if not self.covels_score is None:       comment+="Covels_score="+str(self.covels_score)+" "
    if not self.secisearch_pattern is None: comment+="SECISearch1_pattern="+str(self.secisearch_pattern)+" "
#    comment+='grade='+str(self.grade())+ " energy="+str(self.energy)+' '
#    if not self.type is None:    comment+=" type="+str(self.type)+' '
#    comment+='seq='+nogap(self.alignment.seq_of('t'))+' ss='+self.out_ss(aligned=False)+' '
    return gene.gff(self, tag='secis', program='Seblastian', id='SeB.'+self.id, comment=comment[:-1])

#########           
######## some functions  used for SeB_blast and SeB_exonerate

def p2g_assign_label(self):
    self.label= self.compute_label()
def p2g_blast_evalue(self):
  if isinstance(self, SeB_exonerate): return self.seed.evalue
  else:                               return self.evalue
def p2g_fasta_header(self):  return 'SeB.'+self.id+' label:'+self.label+' class:'+str(self.__class__.__name__)+' '+self.header(no_id=True, no_species=True, no_target=True)
def p2g_one_line_summary(self):  return '-ID:'+self.id.ljust(10)+' '+self.label.ljust(17)+' '+str(self.__class__.__name__).ljust(15)+'   chromosome:'+self.chromosome+' strand:'+self.strand+' '+self.positions_summary(compress=True, max_exons_uncompressed=3)

# all_labels = ['empty','secis_in_cds','secis_not_downstream','pseudo','unaligned', 'lost_sec_candidate','no_sec_candidate','other_ugas','knownSP','newSP']
def p2g_compute_label(self):
    """ function used to compute the label assigned  to these objects"""
    if not self:        return 'empty' 
    if self.overlaps_with(self.secis, phase=False): return 'secis_in_cds'
    if self.distance_from_secis() is None:          return 'secis_not_downstream'
    if 'u' in self.alignment.seq_of('t'):  return 'other_ugas'
    if '*' in self.alignment.seq_of('t') or self.frameshifts():     return 'pseudo'

    if isinstance(self, SeB_exonerate) and 'U' in self.seed.protein():    ## exonerate whose blast had category END will not enter here; also, blasth passing keep_blast will not enter
      ## now filling sec_codons and checking that at least one overlaps with the current prediction.
      sec_codons=[];      some_is_overlapping=False; gaps_target=0
      for pos_in_ali, aa_target  in enumerate(self.seed.alignment.seq_of('t')): #parsing original blast to check if the sec_uga candidate(s) is still aligned or not
        if aa_target=='-':     gaps_target+=1
        elif aa_target=='U':   
          c=self.seed.subseq(1+(pos_in_ali-gaps_target)*3, 3, True)
          if self.overlaps_with(c): some_is_overlapping=True
          sec_codons.append(c)
      if not some_is_overlapping:   return 'unaligned'
      elif  not 'U' in self.alignment.seq_of('t'): return 'lost_sec_candidate' #blast had it, exonerate doesn't

    elif not isinstance(self, SeB_exonerate) or not 'U' in self.seed.protein():  #END category
      if not 'U' in self.alignment.seq_of('t'): return 'no_sec_candidate'

    """
    if isinstance(self, SeB_exonerate) and 'U' in self.seed.protein():    ## exonerate whose blast had category END will not enter here; also, blasth passing keep_blast will not enter
      gaps_target=0
      for p, aa_target in enumerate(self.alignment.seq_of('t')):
        if aa_target=='-': gaps_target+=1
        elif aa_target  in 'U':
          this_codon=self.subseq(1+(p-gaps_target)*3, 3, True)
          if not any([this_codon.overlaps_with(c)  for c in sec_codons] ) : #this is a UGA, but doens't overlap with any of those previously called sec_uga, but is a new sec candidate
            return 'changed_sec_candidate'"""
    if self.category in ['QUE_U', 'ALI_U']: return 'knownSP'
    else: return 'newSP'

def p2g_output_fasta(self):
  """ for outputing"""
  o=">"+self.fasta_header()+'\n'
  o+=fasta(self.protein())
  return o

def p2g_output_cds(self):
  """ for outputing"""
  o=">CDS_of_"+self.fasta_header()+'\n'
  o+=fasta(self.cds())
  return o

def p2g_three_prime(self, length=None):
  """ Returns the sequence downstream of the result. the length provided is the maximum length returned: if not available, a shorter string is returned (can be empty as well in the extreme case)"""
  if length==None: length=opt['three_prime_length']
  d=self.downstream(0, length)
  d.check_boundaries(chromosome_lengths[d.chromosome])
  if d.length():    seq=upper(d.fasta_sequence()[1])
  else: seq=''
  return seq

def p2g_downstream_until_stop(self):
  try: self.downstream_until_stop_seq
  except: self.downstream_until_stop_seq=None
  if self.downstream_until_stop_seq==None:
    seq=replace(self.three_prime(), 'U', 'T')
    index_codon=0; s=''
    while index_codon*3<len(seq):
      codon=seq[index_codon*3:index_codon*3+3]
      s+=codon
      if codon in ['TAA', 'TGA', 'TAG']: index_codon=len(seq) #if found stop codon, do not elongate anymore
      index_codon+=1
    self.downstream_until_stop_seq=s
  return self.downstream_until_stop_seq
    
def p2g_output_secis(self):  return self.secis.fasta()

def p2g_distance_from_secis(self):
  s= self.secis.is_downstream_of(self)
  if s: return s-1
  else: 
  
    return None

def p2g_output_gff(self):
  """ for outputing"""
  #sec positions is automatically added by gff method if there's a sec_positions_list method
  comment='label:'+self.label
  o=gene.gff(self, comment=comment, program='Seblastian', id='SeB.'+self.id)+'\n'
  o+=self.secis.gff()
  return o

def p2g_protein(self):  return nogap(self.alignment.seq_of('t'))

def p2g_sec_pos_in_seq(self):
  """Returns the list of positions 1-based of the U along the protein sequence """
  out=[]
  seq=upper(self.protein())
  f=seq.find('U')
  while f!=-1:
    out.append(f+1)
    f=seq.find('U', f+1)
  return out

def p2g_sec_positions_list(self):
  """ useful to define it to have the sec positions printed in the gff by the gene.gff method"""
  out=[]
  for p in self.sec_pos_in_seq():
    subseqed=self.subseq((p-1)*3+1, 1)
    out.append(subseqed[0][0])
  return out

def p2g_summary(self):
    """ Text representation for the prediction """
    o='_|*|'*8+' | '+str(self.label).center(20)+' | '+ ' secis_index:'+(self.id).ljust(6)+self.__class__.__name__+' category:'+str(self.category)+' '+'|*|_'*8
    o='_'*33+'.'+22*"_"+"."+"_"*(len(o)-57)+'\n'+o+'\n'+'_'*33+'|'+22*"_"+"|"+"_"*(len(o)-57)+'\n\n'
    try:
      o+=' Blastx evalue: '+str(self.blast_evalue())+'\n'
      if self.query_full_title:    o+=' Query name: '+self.query_full_title+'\n'
      else:                        o+=' Query name: '+self.query.chromosome+'\n'    
      o+=' Query range: '+join([str(i) for i in self.query.boundaries()], '-')+'  (total '+str(self.query_prot_length)+')\n\n'
      o+=' Target: '+self.target+'\n'
      o+=' Chromosome: '+self.chromosome+' Strand: '+self.strand+'\n'
      o+=' Positions:\n'
      for exon_index in range(len( self.exons )):
        p_start, p_end= self.exons[exon_index]
        o+=(' -Exon '+str(exon_index+1)).ljust(10)+str(p_start).ljust(12)+' '+str(p_end)+'\n'
      o+='\n'
      o+=' Redox box: '+{True:'YES', False:'no'}[self.has_redox_box()] +'\n'
      o+=self.pretty_alignment()
      o+=' Sequence until next stop: \n'
      s=self.downstream_until_stop()
      if s: 
        translation= ' '+join(list(transl(s)), '  ')+' '
        o+=s+'\n'+translation+'\n'      
      else: 
        o+=' - No sequence available - '
      o+='\n'
      o+=' Distance CDS-SECIS: '+str(self.distance_from_secis())+' ; UGA-SECIS: '+str(self.distance_uga_from_secis())+'\n\n'
      o+=self.secis.summary()
    except:  
      o+='\nERROR while producing this output! skipping...'
    return o

def p2g_distance_uga_from_secis(self):
  """ return the distance (last) secUGA to SECIS"""
  if not self.sec_pos_in_seq(): 
    return None
  dfs=self.distance_from_secis()
  if dfs:  return  dfs + (len(self.protein())-self.sec_pos_in_seq()[-1])*3
  else: return None

def p2g_has_redox_box(self):
  """ returns True if any of the Sec positions has a C or U at 2 aminoacids of distance"""
  seq=self.protein()
  for pos in self.sec_pos_in_seq(): #pos is 1-based
    if pos>=4 and seq[pos-4] =='C': return True  # CxxU
    elif  len(seq)>pos+2 and   seq[pos+2] in 'CU': return True # UxxC or UxxU
  return False
  
p2g_common_methods=[ 'assign_label','compute_label', 'fasta_header', 'output_fasta','output_gff', 'protein', 'sec_pos_in_seq', 'sec_positions_list', 'one_line_summary', 'output_secis', 'summary', 'distance_uga_from_secis', 'has_redox_box', 'distance_from_secis', 'blast_evalue', 'output_cds', 'three_prime', 'downstream_until_stop']
####
class SeB_blast(blasthit):
    """ blast hit subclass adapted for Seblastian. see methods at selenoprofiles_2.py and here below. attributes as in selenoprofiles_2.py plus:
    self.secis=secis class object
    self.category=string, defined in runtime (see code)
    self.blast_query=gene object originally used as query
    """
    def __init__(self, bline=False):
      """ Can be initialized with a simple blaster line, a blaster line with preceded by the category and a tab (in both these cases it is considered that the parsed blast line comes from the blastx of the upstream sequences, and they will be restored), or the output of the abstract line method (no restore needed)"""
      blasthit.__init__(self)
      self.cds_sequence=None
      self.secis=None; self.category=None; self.query_prot_length=None; self.upstream_length=None; self.blast_query=None; self.evalue=None
      self.query_full_title=None
      self.query_full_seq=None

      if bline:
          if len( bline.split('\t') )<=2: self.load_blaster_line(bline, restore_coord=True)
          else:                           self.load_abstract_line(bline)

    def load_blaster_line(self, bline, restore_coord=False, secis_index=False):
      """ Initialized with a blaster parser line (parsed with options -v TLENGTH=1 -v QLENGTH=1) or with the same line preceded by the category and a tab separator"""
      splt= bline.split('\t')
      if len(splt)>1: 
        self.category=splt[0]
        bline=splt[-1]
      else: self.category=None
      try:
        self.load_blaster(bline) 
      except: 
        write('WARNING skipping bad blast hit: '+bline, 1)
        return None

      if not secis_index:      secis_index=int(self.chromosome.split('upstream_of_SECIS_')[1])
      self.id=str(secis_index)
      if not self.secis: self.secis = get_secis(secis_index) ##### loadin .secis if it's not there already (as with load_abstract_line)
      self.query_prot_length = int(bline.split()[11])
      self.upstream_length=     int(bline.split()[12])
      self.blast_query=self.secis.upstream(length_gap, self.upstream_length)
      self.target=target_file
      if restore_coord:      self.restore_absolute_coordinates(self.blast_query)

    def load_blaster(self, line, id=None):
      """ overriding the function in blasthit class to have a reversed blast hit, meaning, we read a blastx hit but we load a tblastn like"""
      splt=line.split()
      self.chromosome=splt[4]
      self.strand='+' 
      start, stop= int(splt[5]), int(splt[6])
      self.add_exon(start, stop)
      #alignmentcd 
      self.alignment=alignment()
      self.alignment.add('q', splt[7])
      self.alignment.add('t', splt[8] ) #not keeping the title here, I already have it in chromosome
      self.evalue=e_v(splt[9])
      if not id is None:      self.id=id
      else:                   self.id= str(uniq_id(self))
      #query
      start, stop= int(splt[1]), int(splt[2])
      self.query=gene(chromosome=splt[0], id=self.id+'_query')
      self.query.strand='+'
      self.query.add_exon(start, stop)

    """
    def abstract_line(self):
      " one line abstract of all information in the object:   id TAB category TAB blasterline TAB secis_abstract_line"
      o=self.id+'\t'+str(self.category)+'\t'
      o+=self.chromosome+' '+str(self.boundaries()[0])+' '+str(self.boundaries()[1])+' + '+self.query.chromosome+' '
      if    self.query.strand=='+':        o+=str(self.query.boundaries()[0])+' '+str(self.query.boundaries()[1])  
      elif  self.query.strand=='-':        o+=str(self.query.boundaries()[1])+' '+str(self.query.boundaries()[0])  
      o+=' '+self.alignment.seq_of('t')+' '+self.alignment.seq_of('q')+' '+str(self.evalue)+' '+str(int( +self.alignment.sequence_identity()*100 ))+' '+str(self.target_length)+' '+str(self.query_length)
      o+='\t'+self.secis.abstract_line()
      return o 

    def load_abstract_line(self, bline):
        "Load from the summary returned by abstract line method "
        self.__init__()
        splt= bline.split('\t')
        self.secis=infernalhit(); self.secis.load_abstract_line(splt[3])
        self.load_blaster_line(splt[2], secis_index=int(splt[0]))
        self.category=splt[1]            
    """

    def set_ugas(self, uga_list):
      """The uga list is a hash of positions in absolute coordinates corresponding to the central nt of the uga codon"""
      try:
        seq=self.alignment.seq_of('t')
        gaps_target=0
        for p_ali, aa_target in enumerate(seq):
          if seq[p_ali]=='-': gaps_target+=1
          elif aa_target=='*':
            this_codon=self.subseq(1+(p_ali-gaps_target)*3, 3)
            is_uga = uga_list.has_key(this_codon[0][0]+1)
            if is_uga: seq=seq[:p_ali]+'u'+seq[p_ali+1:]
        self.alignment.set_sequence('t', seq)
      except:
        print self.query.chromosome
        print self.query
        print self
        raise

    def assign_category(self):
      """Routine to assign blast category. Uses global variable proteins_with_annotated_u  and max_unaligned_end """      
      if proteins_with_annotated_u.has_key(self.query.chromosome):            
        if any(  [ self.alignment.seq_of('t')[p]=='u' and self.alignment.seq_of('q')[p]=='U'   for p in range(self.alignment.length()) ]  ): self.category= 'ALI_U'
        else:                                                                                                                                self.category= 'QUE_U'
      elif any(  [ self.alignment.seq_of('t')[p]=='u' and self.alignment.seq_of('q')[p]=='X'   for p in range(self.alignment.length()) ]  ): self.category= 'ALI_X'
      elif any(  [ self.alignment.seq_of('t')[p]=='u' and self.alignment.seq_of('q')[p]=='C'   for p in range(self.alignment.length()) ]  ): self.category= 'ALI_C'
      elif self.query_prot_length - self.query.boundaries()[1] <= max_unaligned_end:                                                         self.category= 'END'
      else:                                                                                                                                  self.category= 'NOT'
      #now filling U in sequence
      if self.category.startswith('ALI_'):
        aa_category=self.category[-1]
        gaps_target=0
        for p in range(self.alignment.length()):
          aa_query=self.alignment.seq_of('q')[p]; aa_target=self.alignment.seq_of('t')[p]
          if aa_target=='-': gaps_target+=1 #not used actually
          elif aa_target=='u' and aa_query==aa_category:
            self.alignment.set_sequence('t',    self.alignment.seq_of('t')[:p] + 'U' + self.alignment.seq_of('t')[p+1:]  )

    def cds(self, **args): #args is just to pass silencely the useless arguments from aligned_cds function in selenoprofiles
      """ Returns the coding sequence computed with lazy loading"""
      if not self.cds_sequence: 
        self.cds_sequence=upper(self.fasta_sequence()[1])
      return self.cds_sequence


    def get_query_full_title(self):
      """ lazy computing approach to get the full title from the protein database"""
      if self.query_full_title is None: 
        self.fill_full_query_information()
      return self.query_full_title

    def get_query_full_sequence(self):
      """ same for sequence, actually they both use the same subroutine"""
      if self.query_full_seq is None:
        self.fill_full_query_information()
      return self.query_full_seq

    def fill_full_query_information(self):
      """ used by the two previous functions, get_query_full_title and get_query_full_sequence"""
      if self.query_full_seq is None or self.query_full_title is None:
        self.query_full_title, self.query_full_seq = get_sequence_from_protein_db(self.query.chromosome)

    def remove_Xs(self):
      """ It removes Xs from the query and target sequence. The Xs in the query may be coming from the joining p2g procedure (x) or from the blast filtering of query. the target sequence may contain X from the blast filtering(masking). Note that Xs may be present in the original sequences as well, in this case we can't recover them """
      try:
        while 'x' in self.alignment.seq_of('q'):
          full_query_sequence= self.get_query_full_sequence()
          seq=self.alignment.seq_of('q'); seq_no_gap=nogap(seq)
          start_match=seq_no_gap.find('x'); match_lengh=1
          while seq_no_gap[start_match+match_lengh]=='x': match_lengh+=1
          missing_subseq=full_query_sequence[self.query.boundaries()[0]+start_match-1 : self.query.boundaries()[0]+start_match-1+match_lengh]
          start_match_in_ali=self.alignment.position_in_ali('q', start_match+1 ) -1 #0 based now
          gaps_found=0; i=0
          while   i-gaps_found < match_lengh :
            if seq[start_match_in_ali+i]=='-':            gaps_found+=1
            else:  seq=seq[:start_match_in_ali+i]+  missing_subseq[i-gaps_found]   +seq[start_match_in_ali+i+1:]
            i+=1
          self.alignment.set_sequence('q', seq)

        if 'X' in self.alignment.seq_of('q'):
          full_query_sequence= self.get_query_full_sequence()
          start_match_in_ali=self.alignment.seq_of('q').find('X')
          while start_match_in_ali!=-1:
            seq=self.alignment.seq_of('q'); seq_no_gap=nogap(seq)
            start_match=self.alignment.position_in_seq('q', start_match_in_ali+1)-1 #0 based now
            match_lengh=1
            while seq_no_gap[start_match+match_lengh]=='X': match_lengh+=1
            missing_subseq=full_query_sequence[self.query.boundaries()[0]+start_match-1 : self.query.boundaries()[0]+start_match-1+match_lengh]
            gaps_found=0; i=0
            while   i-gaps_found < match_lengh :
              if seq[start_match_in_ali+i]=='-':            gaps_found+=1
              else:  seq=seq[:start_match_in_ali+i]+  missing_subseq[i-gaps_found]   +seq[start_match_in_ali+i+1:]
              i+=1
            self.alignment.set_sequence('q', seq)
            start_match_in_ali=self.alignment.seq_of('q').find('X', start_match_in_ali+i)

        if 'X' in self.alignment.seq_of('t'):
          start_match_in_ali=self.alignment.seq_of('t').find('X')
          while start_match_in_ali!=-1:
            seq=self.alignment.seq_of('t'); seq_no_gap=nogap(seq)
            start_match=self.alignment.position_in_seq('t', start_match_in_ali+1)-1 #0 based now
            match_lengh=1
            while seq_no_gap[start_match+match_lengh]=='X': match_lengh+=1
            gene_subseq=self.subseq(start_match*3+1, match_lengh*3 )
            missing_subseq=    transl(gene_subseq.fasta_sequence()[1])
            gaps_found=0; i=0
            while   i-gaps_found < match_lengh :
              if seq[start_match_in_ali+i]=='-':            gaps_found+=1
              else:  seq=seq[:start_match_in_ali+i]+  missing_subseq[i-gaps_found]   +seq[start_match_in_ali+i+1:]
              i+=1
            self.alignment.set_sequence('t', seq)
            start_match_in_ali=self.alignment.seq_of('t').find('X', start_match_in_ali+i)
      except: 
        printerr('remove_Xs skipping ERROR on this result: '+self.header(), 1)
        


    for method in p2g_common_methods:    exec(method+'= p2g_'+method) 

class SeB_exonerate(exoneratehit):
  """This class is child of the class exoneratehit from selenoprofiles. 
      distance_uga_secis(): returns None if secis is not downstream, the real distance otherwise  (0 if adjacent )
  """
  def __init__(self, exonerateh=None):    
    exoneratehit.__init__(self)
    if not exonerateh is None: self.__dict__=exonerateh.copy(newclass=SeB_exonerate).__dict__
    self.label=None
    self.secis=None
    self.seed=None #blastxhit
    self.cds_sequence=None
    self.query_full_title=None
    self.query_full_seq=None
    self.query_prot_length=None

  def distance_from_secis(self):
    """ returns None if secis is not downstream, the real distance otherwise  (0 if adjacent ) """
    s= self.secis.is_downstream_of(self)
    if s: return s-1
    else: 
      return None

  def set_ugas(self):
    """Same concept as the function for blast hits, but this doesn't require the uga list and is a bit slower"""
    if self:
      seq=self.alignment.seq_of('t')
      gaps_target=0
      for p_ali, aa_target in enumerate(seq):
        aa_query= self.alignment.seq_of('q')[p_ali]
        if aa_target=='-': gaps_target+=1
        elif aa_target=='*':
          #this_codon=self.subseq(1+(p_ali-gaps_target)*3, 3)
          #codon_seq= upper(  this_codon.fasta_sequence()[1]  )
          codon_seq= self.cds()[(p_ali-gaps_target)*3:(p_ali-gaps_target)*3+3]
          if codon_seq=='TGA':           
            if  aa_query in 'CUX': seq=seq[:p_ali]+'U'+seq[p_ali+1:]
            else:                  seq=seq[:p_ali]+'u'+seq[p_ali+1:]
      self.alignment.set_sequence('t', seq)

  for method in p2g_common_methods:    exec(method+'= p2g_'+method)
  __repr__=summary
#######

  """
  def near_pseudos(self):
    " Parse the exonerate file to look for for neighbour exonerate predictions which looks like pseudos, and returns their number 
    "
    pseudo_cont=0
    exonerate_lines=bbash('parse_exonerate.py '+self.ref.file+' -all -seq -l').split('\n')
    for eline in exonerate_lines:
      e=SeB_exonerate()
      e.load(eline, is_string=1)
      if e.is_pseudo():
        pseudo_cont+=1
    return pseudo_cont
    """
#######################################################################################################################################
def transform_secis_ps(infile):
  """ Transforming ps file output by RNAplot to have a nicer image: we rotate and mirror everything. The resulting image still needs to be rotated 180 degrees, which I do with the program convert while producing a png.   I also changed the color of the underlying line"""
  file_h= open(infile, 'r')

  out=''
  line=file_h.readline()
  while line and not line.startswith('/coor'):
    if line=='    dup sequence exch 1 getinterval cshow\n':
      out+='180 rotate\n'+line+'180 rotate\n'    ## rotating letters
    elif  line =='/outlinecolor {0.2 setgray} bind def\n':
      out+='/outlinecolor    {0.75 setgray} bind def\n'
    elif  line =='/paircolor    {0.2 setgray} bind def\n':
      out+='/paircolor    {0.75 setgray} bind def\n'
    else:    out+=line
    line=file_h.readline()
  if not line:  raise Exception, "ERROR format not recognized, couldn't find line /coor"
  out+=line
  line=file_h.readline() #now on first coord line
  coords=[]
  while not line.startswith(']'):
    s=line[1:-2].split()
    x=float(s[0])
    y=float(s[1])
    coords.append([x, y])
    line=file_h.readline()
  min_x = min ( [c[0] for c in coords] )
  max_x = max ( [c[0] for c in coords] )
  mid_x = (min_x+max_x)/2

  def mirror(x, y):    return 2*mid_x-x, y
  for c in coords: #setting mirrored coordinates
    mirrored_c=mirror(c[0], c[1])
    out+= "["+str(round(mirrored_c[0], 3))+' '+str(round(mirrored_c[1], 3))+"]\n"
  while line: #rest of lines (after coordinates)
    out+=line
    line=file_h.readline()

  file_h.close()
  return out

def close_program():
  if 'temp_folder' in globals():
    if 'opt' in globals() and opt['debug']:
      if sys.exc_info()[1]!= None:
        printerr('error preview: '+str(sys.exc_info()[1]), 1)
        raw_input('some error occured. check temp folder '+temp_folder)
    try:     bash('rm -r '+temp_folder)
    except:  pass

  if sys.exc_info()[0]: #an exception was raised. let's print the information using printerr, which puts it in the logfile as well, if any.                                                                                      
    if issubclass(sys.exc_info()[0], SystemExit):      pass
    else:     
      printerr('ERROR '+ traceback.format_exc( sys.exc_info()[2]) , 1)
      sys.exit(8)

  if 'log_file' in globals(): log_file.close()
  try:
    if get_MMlib_var('printed_rchar'):       printerr('\r'+printed_rchar*' ' ) #flushing service msg space
  except:    pass


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except :
    close_program()

