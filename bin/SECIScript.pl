#!/usr/bin/perl
#use CGI qw(:standard);

#NEW CODE FOR SCRIPT VERSION
use threads ;

if ($#ARGV > 10) {
	error_msg() ;
}


use File::Basename;
my $dirname = dirname(__FILE__);

#$patternpath = "/users-d3/mmariotti/software/Seblastian/SSearch" ; #change this to point to the patter files
$temp_folder= "/home/mmariotti/temp"; 				#MM INSERTED CODE
#$bin_folder= "/users-d3/mmariotti/software/Seblastian/SSearch";		#MM INSERTED CODE

### New set up 2021 ## no need to built-in changes anymore!
$patternpath = $dirname . "/SSearch" ;   
$bin_folder=   $patternpath;		


$print_structure = 1;
$draw_image=0;
$totalthreads = 0 ;
$threadnum = 0 ;
$c = '' ; #Default is to not look @ the complement strand
$pat = 'pat_dm2' ; $pattern='default';
%e1=qw(pat_Sep20 -3.7 pat_dm2 -3.7 pat_dm -3.7 pat_s -8.5 pat_c -3.7 pat_g -3.7);
%e2=qw(pat_Sep20 -12.6 pat_dm2 -12.6 pat_dm -11.5 pat_s -15 pat_c -12.6 pat_g -12.6);
$tmpe1 = $e1{$pat} ;
$tmpe2 = $e2{$pat} ;
#print STDERR "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n" ;
#print STDERR "Parameters from command line\n" ;
#print STDERR "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n" ;
for ($i=0; $i <= $#ARGV; $i++) {
#	print "- $ARGV[$i] ($i) -\n" ;
	if (index($ARGV[$i], '-') == 0) { #'cause the i++ thing didn't work
	
		if ($ARGV[$i] eq '-f') { #File name
			$fp = $ARGV[$i+1] ;
			print STDERR "Using sequence file $fp\n" ;
		} elsif ($ARGV[$i] eq '-p') { #Pattern type
			$i = $i + 1 ; #Move to next arg
			$pattern=$ARGV[$i];
			if ($pattern eq 'loose+') { #Loose (canonical and non-canonical)
				$pat = 'pat_c' ;
				print STDERR "Pattern: Loose (canonical and non-canonical)\n" ;
			} elsif ($pattern eq 'loose') { #Loose (canonical)
				$pat = 'pat_Sep20' ;
				print STDERR "Pattern: Loose (canonical)\n" ;
			} elsif ($pattern eq 'strict') { #Strict
				$pat = 'pat_s' ;
				print STDERR "Pattern: Strict\n" ;
			} elsif ($pattern eq 'new') { #New GGA
				$pat = 'pat_g' ;
				print STDERR "Pattern: GGA\n" ;
			} elsif ($pattern eq 'cilia') {
				$pat = 'cilia' ;
				print STDERR "Pattern: ciliates\n" ;
 			} elsif ($pattern eq 'default') {
				print STDERR "Pattern: default\n" ;
			} else {
				die( "Pattern not recognized: $pattern \n" );
			}
	
			$tmpe1 = $e1{$pat} ;
			$tmpe2 = $e2{$pat} ;
			
		} elsif ($ARGV[$i] eq '-complement') {
			$c = '-c' ;
			print STDERR "Complementary strand being searched\n" ;
		} elsif ($ARGV[$i] eq '-temp') {
  		  $temp_folder =  $ARGV[$i+1];
	
	    #print "temp_Folder: *$temp_folder* \n"; 
	
	
		} elsif ($ARGV[$i] eq '-patterns') {
  		  $patternpath =  $ARGV[$i+1];
		} elsif ($ARGV[$i] eq '-bin_folder') {
  		  $bin_folder =  $ARGV[$i+1];
    } elsif ($ARGV[$i] eq '-energy') {
			$i += 1 ; #Jump to first argument
			$tmpe1 = $ARGV[$i] ;
			$i += 1 ;#Jump to second argument
			$tmpe2 = $ARGV[$i] ;
			print STDERR "Custom energy:\n" ;
			print STDERR "    Core: $tmpe1\n" ;
			print STDERR "    Overall: $tmpe2\n" ;
		
		} elsif ($ARGV[$i] eq '-struct') {
			$i += 1 ; #Go to next argument
			print STDERR "Structural filters turned on:\n" ;	
			if (index($ARGV[$i], 'Y') != -1) {
				$Y_flag = 1 ;
				print STDERR "    Y-Filter\n" ;
			}
			if (index($ARGV[$i], 'O') != -1) {
				$O_flag = 1 ;
				print STDERR "    O-Filter\n" ;
			}
			if (index($ARGV[$i], 'B') != -1) {
				$B_flag = 1 ;
				print STDERR "    B-Filter\n" ;
			}
			if (index($ARGV[$i], 'S') != -1) {
				$S_flag = 1 ;
				print STDERR "    S-Filter\n" ;
			}
		} elsif ($ARGV[$i] eq '-image') {
      
      if ($ARGV[$i+1] && ! ($ARGV[$i+1] =~ /^\-/) ){ 
        $i++;
        $draw_image=$ARGV[$i]; 
        }
      else{
        $draw_image=1; 
      }


      
		}
		
	
		elsif ($ARGV[$i] eq '-t') {
			$totalthreads = $ARGV[$i+1] ;
			print STDERR "Using $totalthreads threads\n" ;
		} else {
			error_msg() ;
		}
	}
}



#print STDERR "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-\n\n" ;
$e1 = $tmpe1 ;
$e2 = $tmpe2 ;
#exit ;
#END OF NEW CODE FOR SCRIPT VERSION

#$pat=param('pattern'); #Type of pattern to look for
#$seq=param('sequence'); #Sequence from texarea on the page
#%e1=qw(pat_Sep20 -3.7 pat_dm2 -3.7 pat_dm -3.7 pat_s -8.5 pat_c -3.7);
#%e2=qw(pat_Sep20 -12.6 pat_dm2 -12.6 pat_dm -12.6 pat_s -15 pat_c -12.6);
#$c = param('complstrand') ? '-c' : ''; #search complimentary strand
#$e1=(param('customenergy') && param('e1')) || $e1{$pat}; 
#if (param('customenergy')) { #custom energy parameters from webpage
#	$e1 = param('e1') ;
#} else {
#	$e1 = $e1{$pat} ;
#}
#$e2=(param('customenergy') && param('e2')) || $e2{$pat};
#if (param('customenergy')) {
#	$e2 = param('e2') ;
#} else {
#	$e2 = $e2{$pat} ;
#}

#$Y_flag=param('Y_filter');
#$O_flag=param('O_filter');
#$B_flag=param('B_filter');
#$S_flag=param('S_filter');


#if ($file = param('seqfile')) { #In case a file is to be read
#	$seq=join "", <$file> 
#} 

#if (! $seq) { abnormal("did not get seq") } #Couln't read file

###--sequence reformatting--
#$seq=~s/\r/\n/g; #chanching windows curret returns to UNIX new lines
#@seqarr=split /\n/, $seq;
#for $line (@seqarr) {
#	if ($line=~/^>/) { #removing digits and spaces
#		$line=~s/ |\,/_/g
#	} else {
#		$line=~s/\d| //g
#	} 
#}
#$seq=join "\n", @seqarr;
#unless ($seq =~ /^>/) {
#	$seq=">".( param('name') || "unnamed")."\n".$seq
#}

###--running scan_for_matches with choosen pattern
#open (SEQ, ">../tmppng/$$"."_seq");
#print SEQ $seq;
#close (SEQ);
#print "$fp\n" ;
if (! $fp){die("ERROR you must specify a fasta sequence file with option -f ")}

if ($draw_image == 1){ $image_prefix=$fp. ".secis"}
elsif ($draw_image){   $image_prefix=$draw_image}

print STDERR "Script starting:\n" ;
print STDERR "Running PatScan on sequence file\n" ;
$comm="$bin_folder/scan_1GB $c $patternpath/$pat < $fp  > $temp_folder/$$"."_hits 2> $temp_folder/scanerrlog_SECISearch"; #MM INSERTED CODE
#print STDERR $comm . "\n";
`$comm`;
print STDERR "PatScan run done, moving on\nChecking energy on hits from PatScan  " ;
###--running show_hits

my $hits_file="$temp_folder/$$"."_hits";
open (HITS, $hits_file); #MM INSERTED CODE
$hitnumber=0;
while (<HITS>) {

	if (($hitnumber%4) == 0) {
		print STDERR "\b\\" ;
	} elsif (($hitnumber%4) == 1) {
		print STDERR "\b|" ;
	} elsif (($hitnumber%4) == 2) {
		print STDERR "\b/" ;
	} else {
		print STDERR "\b-" ;
	}

	@param=();
	my $i = 0 ;
	my @threads ;

#	($g1, $name, $start, $stop) = split />|\:\[|\,|\]/;
	/>(.+)\:\[(\d+),(\d+)\]/ ; 				 #MM INSERTED CODE
	($g1, $name, $start, $stop) = ("", $1, $2, $3) ; 	 #MM INSERTED CODE

	$_=<HITS>;
	tr/Tt/Uu/;
	@a=split;
	s/\w/./g for @b=@a;
	$b[3]='x(((('; $b[5]='xx'; 
	if (length ($a[8]) == 4)  { 
		$b[8]='))))' 
	} elsif (length ($a[8]) == 5)  { 
		$b[8]='.))))' 
	} else { 
		#die "\nReal $a[8] - a8 length problem\nIn $name\n" 
	}

	if (($upstemen = energy(3..8)) < $e1 && ($fullstructen = energy(0..11)) < $e2 && structure_OK()) { 
		push @param, $currentseq;
		push @param, $currentstruct;
		push @param, 15;
		push @param, $image_prefix.($hitnumber+1).".png";
		push @param, 170;
		push @param, 250;
		push @param, length (join "", @a[0..2])+1;
		push @param, length (join "", @a[0..2])+4;
		push @param, (length (join "", @a[0..2]))..((length (join "", @a[0..2]))+4);
		push @param, (length (join "", @a[0..4]))..((length (join "", @a[0..4]))+1);
		
		if (length ($a[8]) == 4) { 
			push @param, (length (join "", @a[0..7]))..((length (join "", @a[0..7]))+3) 
		} else { 
			push @param, (length (join "", @a[0..7])+1)..((length (join "", @a[0..7]))+4) 
		}
#Threaded code

		if ($totalthreads > 1) {
		
			@thread = threads->list ;
#		print $#thread + 1 ."\n";
			$threadnum = $#thread +1;
	
			if ($threadnum < $totalthreads) { #We haven't reached the thread quota
				$thr = new threads \&imager2, @param ;
#				print " New Thread created (".$thr->tid.")" ;
			} else {
				if (($hitnumber-$totalthreads) == 0) {
					$thr2 = $thread[0] ;	
#					print "ZERO " ;
				} else {
					$thr2 = $thread[1] ;
#					print "NOT ZERO " ;
				}
				#$thr2 = $thread[$hitnumber-$totalthreads] ;
				$thr = new threads \&imager2, @param ;
#				print "New Thread created (".$thr->tid.")" ;
				$thr2 ->join ;
#				print "Joined thread ".($hitnumber-$totalthreads+1)."\n" ;
			}
#Threaded code
		} else {	
			#print "Calling imager\n" ;
			imager2(@param);
		}
		
		push @texta, [ @a ];
		push @name, $name;
		push @start, $start;
		push @stop, $stop;
		push @png, "\"".$$.$hitnumber.".png\"";
		push @comment, $start>$stop ? "(complement)" : "";
		push @upstemen, $upstemen;
		push @fullstructen, $fullstructen;
		push @structures, $currentstruct;

		$hitnumber++;

	}
}

if ($totalthreads > 1) {
	@thread = threads->list ;
#	print "Exiting the last $totalthreads threads of the bunch\n" ;
	$thr2 = $thread[0] ;
	$thr2->join ;
	$thr2 = $thread[1] ;
	$thr2->join ;
}
print STDERR "\nDone!\n" ;
close(HITS); 
unlink( $hits_file);
results_output();

########################################################################################################
########################################################################################################
#usage energy({list of numbers of PatScan pattern fields (starting from 0)})
#example: energy(0..11) calculates energy for the whole structure in case of classical pattern
#example: energy(3..8) calculates energy for the core SECIS ATGAN..NGAN in case of classical pattern
sub energy {
	open(TMP, ">$temp_folder/$$"."_en");
	print TMP  ">$$\n", @a[@_], "\n", @b[@_];
	close(TMP);
	($g1, $currentseq, $currentstruct, $energy)=split / \(|\n/, `$bin_folder/RNAfold13 -C < $temp_folder/$$'_en'`; #MM INSERTED CODE
	chomp $energy; chop $energy;
	unlink ("$temp_folder/$$"."_en");
	unlink ("$$"."_ss.ps");
	return $energy;
}
##########################################################################################################
sub structure_OK {
	my $fout=1;

	my $substruct=substr($currentstruct, length( join '', @a[0..3]),  length( join '', @a[4..7], (length ($a[8]) == 4) ? () : substr($a[8], 0, 1)   )) ;

#Y-filter detects Y-shaped SECISes
	if (($substruct =~ /\)\.*\(/) && $Y_flag) { 
		$fout=0 
	}

#O-filter detects SECISes with more than 2 umparied nucleotides in a row on any strand in the first 7nt from the quartet
	my $porog=7;
	if ( (substr($substruct, 0, $porog) =~ /\.{3,}/ ||  substr($substruct, -$porog)=~ /\.{3,}/ ) && $O_flag) { 
		$fout=0 
	}

#B-filter detects SECISes that have more than 2 unpaired nt more on 5' side that on 3' side (visually bended to the left)
	($substruct =~ /(.*)(\(\.*\))(.*)/);
	my $up=$1; 
	my $down=$3;

	if (((($up =~ s/\./\./g)-($down =~ s/\./\./g)) < (-2)) && $B_flag) { 
		$fout=0 
	}

#S-filter detects structures with less that 8 pairs
	if ((( $substruct =~ tr/\)/\)/ ) < 8 ) && $S_flag ) { 
		$fout=0 
	}

	return $fout;
}

######################Imager2 block start################################################################################
sub imager2 {
	@bold=();
	(my $seq, my $struct, $scale, $outfilename, $offx, $offy, $posbase1, $posbase2, @bold)=@_;

	@x=();
	@xnew=();
	@ynew=();
	@y=();
	@seq=();
	@xf=();
	@yf=();

	$pi=3.1415926;
	$a=$scale;
	$fontsize=$a;
	$fonttype='Courier';
	$seq=' '.$seq.' ';
	$structmod="(".$struct.")";

	@seq=split '', $seq;

	$x[0]=0;
	$y[0]=0;


	$x[$#seq]=$x[0]+$a;
	$y[$#seq]=$y[0];


	@pairof=pairof_from_brackets($structmod);

	@parts=(0);
	while (@parts) { 
		drawpart(pop @parts) 
	};
	
	pop @x; 
	pop @y; 
	pop @seq; 
	shift @x; 
	shift @y; 
	shift @seq;
	
	@pairof=pairof_from_brackets($struct);
	position_picture($posbase1, $posbase2, $offx, $offy);

  
#UNCOMMENT this for pics
  if ($draw_image){	psoutput(); }

##########################################################################################
sub pairof_from_brackets {
	my @pairof;
	my @struct=split "", $_[0];
	for my $i (0..$#struct) {
		if ($struct[$i] eq '.') { 
			$pairof[$i]=-1 
		} elsif ($struct[$i] eq '(') {
			push @left, $i
		} elsif ($struct[$i] eq ')') {
			my $k=pop @left; 
			$pairof[$k]=$i; 
			$pairof[$i]=$k
		} else {
			die "error: unrecognized symbol \"$struct[$i]\" in structure"
		}
	}
	return @pairof
}
##########################################################################################
sub drawpart {
	$i=$_[0];
	$j=$pairof[$i];
	$k=1;
	@circle=();
	@xlocal=();
	@ylocal=();

	$xlocal[$i]=0;
	$ylocal[$i]=0;

	while ($pairof[$i+$k]==($j-$k)) {
		$xlocal[$i+$k]=0;
		$ylocal[$i+$k]=$a*$k;
		($x[$i+$k], $y[$i+$k])=globalized($xlocal[$i+$k],$ylocal[$i+$k]);
 
		$xlocal[$j-$k]=$a;
		$ylocal[$j-$k]=$a*$k;
		($x[$j-$k], $y[$j-$k])=globalized($xlocal[$j-$k],$ylocal[$j-$k]);
		$k++; 
	}

	for ($n=$i+$k; $n<=$j-$k; $n++) {
    	if ($pairof[$n]!=(-1)) { 
        	push @parts, $n;
        	push @circle, $n;
        	$n=$pairof[$n];
    	}
    	push @circle, $n;
	}

	$alpha=2*$pi/(2+@circle);
	$r=$a/(2*sin($alpha/2));
	$xcenter=$xlocal[$i+$k-1]+$r*sin($alpha/2);
	$ycenter=$ylocal[$i+$k-1]+$r*cos($alpha/2);

	for $ci (0..$#circle) {
		$beta=(3/2)*($pi-$alpha)-$alpha*$ci;
		$xlocal[$circle[$ci]]=$xcenter+$r*cos($beta);
		$ylocal[$circle[$ci]]=$ycenter+$r*sin($beta);
		($x[$circle[$ci]], $y[$circle[$ci]])=globalized($xlocal[$circle[$ci]], $ylocal[$circle[$ci]])
	}
##########################
	sub globalized {
		my $xl=$_[0];        
		my $yl=$_[1];
		$rotangle=atan2( $y[$j]-$y[$i], $x[$j]-$x[$i] );
		my $xg=sqrt( $xl**2+$yl**2 )*cos( atan2($yl, $xl)+$rotangle)+$x[$i];
		my $yg=sqrt( $xl**2+$yl**2 )*sin( atan2($yl, $xl)+$rotangle)+$y[$i];
		return ($xg, $yg);
	}
#########################
}
#################################################################################################
sub position_picture {
	my ($base1, $base2, $xbase1offset, $ybase1offset) = @_;
	my $alpha=3.1415926/2-atan2 ($y[$base2]-$y[$base1], $x[$base2]-$x[$base1]) ;

	for $j (0..$#x) {
		$xnew[$j]=$x[$j]-$x[$base1]; 
		$ynew[$j]=$y[$j]-$y[$base1];
		($xnew[$j], $ynew[$j])= (sqrt($xnew[$j]**2+$ynew[$j]**2)*cos( atan2($ynew[$j],$xnew[$j])+$alpha),sqrt($xnew[$j]**2+$ynew[$j]**2)*sin( atan2($ynew[$j],$xnew[$j])+$alpha));
		$xnew[$j]+=$xbase1offset;
		$ynew[$j]+=$ybase1offset;
	}

	@x=@xnew;
	@y=@ynew;

}

#################################################################################################
sub psoutput {
	open(OUTF, "| tee seciscript_log |  /usr/bin/gs -q -dNOPAUSE -sDEVICE=png256 -g400x600 -dBATCH -sOutputFile=$outfilename -");
	show_seq();
	show_connections();
	show_bounds();
	print OUTF "stroke\n showpage\n";
}
############################
sub show_seq {

	@nonbold=(0..$#x);
	$k=0;

	print OUTF "newpath\n /$fonttype"."Bold findfont $fontsize scalefont setfont\n";
	for $index (@bold) {
		$xf=$x[$index]-0.6*$fontsize/2;
		$yf=$y[$index]-0.6*$fontsize/2;
		print OUTF "$xf $yf moveto\n ($seq[$index]) show\n";

		splice @nonbold, $index-$k, 1; 
		$k++;
	}

	print OUTF "newpath\n /$fonttype findfont $fontsize scalefont setfont\n";
	for $index (@nonbold) {
		$xf=$x[$index]-0.6*$fontsize/2;
		$yf=$y[$index]-0.6*$fontsize/2;

		if (! defined $seq[$index]) { 
			die "@seq\n@x\n"
		}
		print OUTF "$xf $yf moveto\n ($seq[$index]) show\n";
	}
}
#############################
sub show_connections {
	for my $i (0..($#x-1)) {
		$x1=$x[$i]; 
		$y1=$y[$i];
		$x2=$x[$i+1];
		$y2=$y[$i+1];
		line_bound('black', 0.25, 0.04); 
	}
}
#############################
sub show_bounds {
	$gam=0.2;

	for my $i (0..$#x) {
        if ($i<$pairof[$i]) {
             $curpair=$seq[$i].$seq[$pairof[$i]];
             $x1=$x[$i]; $y1=$y[$i]; $x2=$x[$pairof[$i]]; $y2=$y[$pairof[$i]];

			if (($curpair eq 'GU') || ($curpair eq 'UG')) {
				 circle_bound('lightblue', $gam)
			} elsif (($curpair eq 'AU') || ($curpair eq 'UA')) {
					circle_bound('mediumblue', $gam)
			} elsif  (($curpair eq 'GC') || ($curpair eq 'CG')){
				circle_bound('red', $gam)
			} else {
				diamond_bound('burgundy', $gam, 0.3)
			}
       }
	}
}
###########################
sub line_bound {
	my @color=rgbcolor($_[0]);
	my $gamma=$_[1];
	my $teta=$a*$_[2]; 

	$xstart=$x1+($x2-$x1)*(1-$gamma)/2; 
	$ystart=$y1+($y2-$y1)*(1-$gamma)/2;
	$xstop=$x2-($x2-$x1)*(1-$gamma)/2;
	$ystop=$y2-($y2-$y1)*(1-$gamma)/2;
	print OUTF "newpath\n @color setrgbcolor\n $teta setlinewidth\n $xstart $ystart moveto\n $xstop $ystop lineto\n stroke\n";
}
###########################
sub circle_bound {
	my @color=rgbcolor($_[0]);
	my $gamma=$_[1];
	my $xc=($x1+$x2)/2;
	my $yc=($y1+$y2)/2;
	my $rc=$a*$gamma/2;
	print OUTF "newpath\n @color setrgbcolor\n $xc $yc $rc 0 360 arc\n fill\n stroke\n";
}
###########################
sub diamond_bound {
	my @color=rgbcolor($_[0]);
	my $gamma=$a*$_[1]/2;
	my $teta=$a*$_[2]/2; 
	my $rotangle=atan2 ($y2-$y1, $x2-$x1);
	my $xc=($x1+$x2)/2;
	my $yc=($y1+$y2)/2;
	my $x_1=(-1)*sin($rotangle)*$teta+$xc;
	my $y_1=cos($rotangle)*$teta+$yc;
	my $x_2=cos($rotangle)*$gamma+$xc;
	my $y_2=sin($rotangle)*$gamma+$yc;
	my $x_3=sin($rotangle)*$teta+$xc;
	my $y_3=(-1)*cos($rotangle)*$teta+$yc;
	my $x_4=(-1)*cos($rotangle)*$gamma+$xc;
	my $y_4=(-1)*sin($rotangle)*$gamma+$yc;
	print OUTF "newpath\n @color setrgbcolor\n $x_1 $y_1 moveto $x_2 $y_2 lineto $x_3 $y_3 lineto $x_4 $y_4 lineto closepath\nfill\n stroke\n";
}
###########################
sub rgbcolor {
	my @rgb;
	if ($_[0] eq 'lightblue') { 
		@rgb=qw(0.7 0.7 0.7) 
	} elsif ($_[0] eq 'mediumblue') { 
		@rgb=qw(0 0 0.9) 
	} elsif ($_[0] eq 'darkblue') { 
		@rgb=qw(0 0 0.2) 
	} elsif ($_[0] eq 'burgundy')  { 
		@rgb=qw(0.5 0.1 0.1) 
	} elsif ($_[0] eq 'black')  { 
		@rgb=qw(0.5 0.5 0.5) 
	} elsif ($_[0] eq 'red')  { 
		@rgb=qw(0.9 0 0) 
	}

	return @rgb;
	}
}
###############Imager2 block stop##########################################################################################


####################################################################################################
sub results_output {
#	print header; 
#	print "<HTML><HEAD><TITLE>SECISearch</TITLE></HEAD><BODY BGCOLOR=CCCCCC>\n";
	$totalhits=$#name+1;
	if (@name) { 
		print "#Pattern: $pattern ; Total hits: $totalhits\n";
#		print "<TABLE BORDER=1>";

		for (0..$#name) {
#			print "-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*\n" ;
			#$name[$_]=~s/_/ /g;
			$hn=$_+1;
			if ($start[$_] < $stop[$_]) {
				print ">$hn $name[$_]: $start[$_] $stop[$_] $comment[$_]\n" ;
			} else {
				print ">$hn $name[$_]: $stop[$_] $start[$_] $comment[$_]\n" ;
			}
			if ($pat eq 'pat_c') {
#				print "upstemen:$upstemen[$_] fullstructen: $fullstructen[$_]\n" ;
			}
      $seq=join(' ', @{$texta[$_]}[0..2]). " ". $texta[$_][3]." ". $texta[$_][4] . " ". $texta[$_][5]. " ". join(' ', @{$texta[$_]}[6..7]) . " " . $texta[$_][8] . " " . join( ' ', @{$texta[$_]}[9..11]) ; 
			print $seq . "\n"; 

		  if ($print_structure){ 
		    %white_space_pos=();
		    $i=0;
		    $i= index($seq, ' ', $i);
		    while ( $i != -1 ) {
		      $white_space_pos{$i}=1;
    		  $i= index($seq, " " , $i+1)
		    }
		
    		$sstructure=''; $spaces=0;
		    for ($i=0; $i< length( $structures[$_]); $i++  ){

		      if ($white_space_pos{$i+$spaces}){
		        $sstructure.=' ';
		        $spaces+=1;
		        }
          $sstructure.=substr($structures[$_],$i, 1);		    
		      }
		
		  #### note the lengths of the parts printed. print ss in consequence
		  print  $sstructure. "\n"; 		
		   }

		}
		
#	print "</TABLE>";
	} else { 
		print "#Pattern: $pattern ; Total hits: 0 ; Sorry, no SECIS elements were found! \n";
#		if ($pat ne "pat_Sep20") { 
#			print "You may want to try the other pattern" 
#		} 
	}
exit ;
}
####################################################################################################
####################################################################################################

#print end_html;


######################################
sub abnormal {
#	print header; 
	print "<HTML><HEAD><TITLE>SECISearch</TITLE></HEAD><BODY BGCOLOR=CCCCCC>";
	print $_[0];
#	print end_html;
	exit;
}

sub error_msg {
	print "Usage: SECIScript -f <FILE> -struct <YOBS> [-p <PATTERN>] [-complement] [-energy <CORE> <OVERALL>] [-t <NUMBER>]\n" ;
	print "		-f <FILE> is the FastA file containing the sequence of interest\n" ;
	print "		-struct [YOBS] are set when Y, O, B and S filters are going to be used.\n" ;
	print "[OPTIONAL parameters]\n" ;
	print "		-p <PATTERN> can be: \n" ;
	print "			'loose+' (Loose cannonical and non-canonical)\n" ;
	print "			'loose' (Loose canonical)\n" ;
	print "			'strict' (Stric)\n" ;
	print "		-complement is used the compliment of the sequence should also be searched\n" ;
	print "		-energy <CORE> <OVERALL> is used when custom energy paramenters are going to be used.\n" ;
	print "			<CORE> is the value of the energy of the core structure (default = -5)\n" ;
	print "			<OVERALL> is the value of the energy of the overall structure (default = -11)\n" ;
	print "		-t <NUMBER> is the number of threads that should be spawned\n" ;
	print "   -image   to produce image as output when one or more SECISes are found\n";
	print "   -temp         provide temporary folder for intermediate files\n";
	print "   -patterns     provide folder where patterns are searched for\n";
	print "   -bin_folder   provide folder where executables (SCAN1Gb, RNAfold) are searched for\n";
	
	
	exit ;
	
}
