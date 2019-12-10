=head1 LICENSE
Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT
 Ensembl <http://www.ensembl.org/info/about/contact/index.html>

=cut

=head1 NAME
 five_prime_UTR_annotator
=head1 SYNOPSIS
 mv five_prime_UTR_annotator.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin five_prime_UTR_annotator, uORF_starts_ends_GRCh37_PUBLIC.txt
=head1 DESCRIPTION
 A VEP plugin that annotates the effect of 5' UTR variant especially for variant creating start codon uAUG and disrupting stop codon
 Please cite Whiffin et al. Characterising the loss-of-function impact of 5' untranslated region variants in whole genome sequence data from 15,708 individuals. bioRxiv (2019)
=cut


package five_prime_UTR_annotator;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

    sub feature_types {
        return ['Transcript'];
    }

	sub new {
	my $class = shift;
	my $self = $class->SUPER::new(@_);
	my $file = $self->params->[0];

	if(!$file) {
    my $plugin_dir = $INC{'five_prime_UTR_annotator.pm'};
    $plugin_dir =~ s/five_prime_UTR_annotator\.pm//i;
    $file = $plugin_dir.'/uORF_starts_ends_GRCh37_PUBLIC.txt';
  }

  die("ERROR: SORF file $file not found\n") unless $file && -e $file;

	open my $fh, "<",  $file;
  	my %stop_evidence;

  while(<$fh>) {
    chomp;
    my ($chr, $pos, $gene, $strand, $type, $stop_pos) = split /\t/;
	my $from="";
    my $to="";
      if ($strand eq 'forward')
       {
          $from = $stop_pos - 2;
           $to = $stop_pos;     }
      elsif ($strand eq 'reverse')
      {
          $from = $stop_pos;
          $to = $stop_pos + 2;
      }

      for (my $n = $from; $n<=$to; $n++)
     {
          my $key = $chr.":".$n; # chr has 'chr' proceeding
          $stop_evidence{$key} = 1;
     }
  }

  close $fh;

  $self->{stop_evidence} = \%stop_evidence;

  return $self;
}

    sub get_header_info {

	$self->{_header_info} = {
	#If the variant creates an ATG in the five prime UTR sequence
        uAUG_KozakContext => "The surrounding Kozak sequence of the uAUG",
        uAUG_KozakStrength => "Strength of the surrounding Kozak consensus of the uAUG",
        uAUG_DistanceToCDS => 'The uAUG distance upstream of the main ORF coding sequence',
        uAUG_FrameWithCDS => 'Frame with respect to the main ORF coding sequence',
        uAUG_InframeStop => 'Whether there is an inframe stop codon with respect to the uAUG within the 5 prime UTR',
        uAUG_DistanceToInframeStop => 'The uAUG distance to the inframe stop codon within the 5 prime UTR',
        uAUG_DistanceFrmCap => 'Distance of uAUG from 5 prime mRNA cap',
        #If the variant disrupts the stop codon of an existing uORF
        uSTOP_AltStop => 'Whether there is another stop codon in the uORF',
        uSTOP_AltStopDistance => 'The distance of alternative stop codon to the lost stop codon',
        uSTOP_AltStopDistanceToCDS => 'The distance of alternative stop codon to CDS',
        uSTOP_FrameWithCDS => 'Frame with respect to the main ORF coding sequence ',
        uSTOP_KozakContext => 'The surrounding Kozak sequence of the uAUG in the existing uORF',
        uSTOP_KozakStrength => 'Strength of the surrounding Kozak consensus of the uAUG in the existing uORF',
        uSTOP_Evidence => 'Whether there is prior evidence of translation documented in sorfs.org',
        #If either is true
        existing_uORFs => 'The number of existing uORFs already within the 5 prime UTR',
        existing_oORFs => 'The number of existing oORFs already within the 5 prime UTR',
        existing_inframeORFs => 'The number of existing inframe ORFs already within the 5 prime UTR',
        };
		return $self->{_header_info};
    }

    sub run {
        my ($self, $tva) = @_;

  #only annotate the effect if the variant is (1)5_prime_UTR_variant
	return {} unless grep {$_->SO_term eq '5_prime_UTR_variant'}  @{$tva->get_all_OverlapConsequences};

	#retrieve the variant info
	my $vf = $tva->variation_feature;
	my $chr = ($vf->{chr}||$vf->seq_region_name);
	my $pos = ($vf->{start}||$vf->seq_region_start);
	my $ref = $tva->base_variation_feature_overlap->get_reference_VariationFeatureOverlapAllele->variation_feature_seq;
	my $alt = $tva->variation_feature_seq;
	my %variant = (
	"chr" => $chr,
	"pos" => $pos,
	"ref" => $ref,
	"alt" => $alt,
	);
  #only annotate if the variant is SNV
  return {} unless (length($ref)==1 && length($alt)==1);

  #retrieve the UTR info: transcript id, strand, five prime UTR sequence, start and end genomic coordinates.
	my $tv = $tva->base_variation_feature_overlap;
	my $t = $tv->transcript;
	my $transcript_id = (defined $t? $t->stable_id: undef);

	#retrieve the gene symbol of the transcript
	my $symbol = $t->{_gene_symbol} || $t->{_gene_hgnc};
	#retrieve the strand of the transcript
	my $tr_strand = $t->strand + 0;
	#retrieve the five prime utr sequence
	my $five_prime_seq = $t->five_prime_utr->seq();

	#Type: five_prime_feature - Bio::EnsEMBL::Feature
	my $UTRs = $t->get_all_five_prime_UTRs();
	my @five_utr_starts;
	my @five_utr_ends;
	foreach  my $utr (@$UTRs){
		my $start = $utr->start();
		my $end = $utr->end();
		push(@five_utr_starts, $start);
		push(@five_utr_ends,$end);
	}

	my @sorted_starts = sort {$a <=> $b} @five_utr_starts;
	my @sorted_ends = sort {$a <=> $b} @five_utr_ends;

	my %UTR_info = (
	"gene" => $symbol,
	"start" => \@sorted_starts,
	"end" => \@sorted_ends,
	"seq" => $five_prime_seq,
	"strand" => $tr_strand,
	);

	my $creatingATG_effect = $self->creating_ATG(\%variant,\%UTR_info);
  	my $removingstop_effect = $self->removing_stop(\%variant,\%UTR_info);
  	my $existing_uORF = $self->count_number_ATG($five_prime_seq);
	my $utr_effect ={%$creatingATG_effect, %$removingstop_effect, %$existing_uORF};
	return $utr_effect? $utr_effect: {};
 }

sub removing_stop{
    #Description: annotate if a five_prime_UTR_varint removes a stop codon
    #Returntype: hashref

    #TODO: check whether it's SNV

    my ($self, $variant_info, $UTR_info) = @_;
    my %flip;
    $flip{'A'}='T';
    $flip{'C'}='G';
    $flip{'G'}='C';
    $flip{'T'}='A';

    my %kozak_strength;
    $kozak_strength{1}='Weak';
    $kozak_strength{2}='Moderate';
    $kozak_strength{3}='Strong';

    my $chr = $variant_info->{chr};
    my $pos = $variant_info->{pos};
    my $ref = $variant_info->{ref};
    my $alt = $variant_info->{alt};

    my $gene = $UTR_info->{gene};
    my @sorted_starts = @{$UTR_info->{start}};
    my @sorted_ends = @{$UTR_info->{end}};
    my $num_exons = @sorted_starts;
    my @sequence = split //, $UTR_info->{seq};
    my $length = @sequence;
    my $strand = $UTR_info->{strand};

    my %existing_uORF = %{$self->existing_uORF(\@sequence)};

    #return annotators
    my $uSTOP_AltStop = "";
    my $uSTOP_AltStopDistance = "";
    my $uSTOP_AltStopDistanceToCDS = "";
    my $uSTOP_FrameWithStart = "";
    my $uSTOP_KozakContext = "";
    my $uSTOP_KozakStrength = "";
    my $current_kozak = "";
    my $current_kozak_strength = "";
    my $uSTOP_evidence = "";

    #indicate whether the variant disrupts a stop codon
    my $flag=0;

    #the relative position of input variant in the UTR sequence
    my $n;
    my %chr_position;
    my $utr_position = 0;

    my $result={};

    if ($strand == 1){

        #create a map from chromosome position to UTR position
        for (my $m=0; $m<$num_exons; $m++)
        {
            for (my $p=$sorted_starts[$m]; $p<=$sorted_ends[$m]; $p++)
            {
                $chr_position{$p}=$utr_position;
                $utr_position++;
            }
        }
        #finding the relative pos of variant in the UTR sequence
        $n = $chr_position{$pos};

        my @start = sort {$a <=> $b} (keys %existing_uORF);
        my $start_pos;


        for ($i=0;$i<=@start;$i++){
            $start_pos = $start[$i];
            my @stops = sort {$a <=> $b} @{$existing_uORF{$start_pos}};
            my $end_pos=$stops[0];

            my $base_pos = $n-$end_pos;
            if (($base_pos>2) || ($base_pos<0)){next;};

            my $stop_codon= $sequence[$end_pos].$sequence[$end_pos+1].$sequence[$end_pos+2];
            #TAA->TAG or TAA->TGA is still a stop codon
            if (($stop_codon eq 'TAA') && ($base_pos) && ($alt eq 'G')){next;}
            #TGA->TAA
            if (($stop_codon eq 'TGA') && ($base_pos == 1) && ($alt eq 'A')){next;}
            #TAG->TAA
            if (($stop_codon eq 'TAG') && ($base_pos == 2) && ($alt eq 'A')){next;}
            #if the variant indeed disrupts a stop codon, stop the loop
            $flag = 1;

            if($flag){

                #getting the Kozak context and Kozak strength of the start codon

  				if ((($start_pos-3)>=0)&&($sequence[($start_pos+3)])){
  					$current_kozak = $sequence[($start_pos-3)].$sequence[($start_pos-2)].$sequence[$start_pos-1]."ATG".$sequence[$start_pos+3];
  				}
  				else{
  					$current_kozak = '-';
  				}

                if ($current_kozak !~ /-/){
                    my @split_kozak = split //, $current_kozak;
                    $current_kozak_strength = 1;
                    if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
                        $current_kozak_strength = 3;
                    }
                    elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
                        $current_kozak_strength = 2;
                    }
                }
                if($current_kozak_strength>$uSTOP_KozakStrength){
                    $uSTOP_KozakStrength = $current_kozak_strength;
                    $uSTOP_KozakContext = $current_kozak;
                }

                #if there is an alternative stop codon in the uORF
                if (@stops>1){
                $uSTOP_AltStop = "True";
                $uSTOP_AltStopDistance = $stops[1]-$stops[0];
                $uSTOP_AltStopDistanceToCDS = $length-$stops[1];
                } #if there is no alternative stop codon
                else{
                    $uSTOP_AltStop = "False";
                    $uSTOP_AltStopDistance = "-";
                    $uSTOP_AltStopDistanceToCDS = "-";
                }
            	  if (($length-$start_pos) % 3){
                	$uSTOP_FrameWithStart = "outOfFrame";
                    }
                else{
                    $uSTOP_FrameWithStart = "inFrame";
                    }
            	#find evidence from sorf
            	my $query = "chr".$chr.":".$pos;
            	$uSTOP_evidence=$self->{stop_evidence}->{$query}?"True":"False";

            }

        }
        }

    if ($strand == -1){
        for (my $m=$num_exons-1; $m>=0; $m--)
        {
            for (my $p=$sorted_ends[$m]; $p>=$sorted_starts[$m]; $p--)
            {
                $chr_position{$p}=$utr_position;
                $utr_position++;
            }
        }
        #finding the relative pos of variant in the sequence (5' to 3' correct strand)
        $n = $chr_position{$pos};

        my @start = sort {$a <=> $b} (keys %existing_uORF);
        my $start_pos;

        for ($i=0;$i<=@start;$i++){
            $start_pos = $start[$i];
            my @stops = sort {$a <=> $b} @{$existing_uORF{$start_pos}};
            my $end_pos=$stops[0];

            my $base_pos = $n-$end_pos;
            if (($base_pos>2) || ($base_pos<0)){next;};
            my $stop_codon=$sequence[$n].$sequence[$n+1].$sequence[$n+2];
            #TAA->TAG or TAA->TGA is still a stop codon
            if (($stop_codon eq 'TAA') && ($base_pos) && ($alt eq 'C')){next;}
            #TGA->TAA
            if (($stop_codon eq 'TGA') && ($base_pos == 1) && ($alt eq 'T')){next;}
            #TAG->TAA
            if (($stop_codon eq 'TAG') && ($base_pos == 2) && ($alt eq 'T')){next;}
            #if the variant indeed disrupts a stop codon, stop the loop
            $flag = 1;

            if($flag){

                #getting the Kozak context and Kozak strength of the start codon

  			if ((($start_pos-3)>=0)&&($sequence[($start_pos+3)])){
  				$current_kozak = $sequence[($start_pos-3)].$sequence[($start_pos-2)].$sequence[$start_pos-1]."ATG".$sequence[$start_pos+3];
  			}
  			else{
  				$current_kozak = '-';
  			}
		     if ($current_kozak !~ /-/){
                    my @split_kozak = split //, $current_kozak;
                    $current_kozak_strength = 1;
                    if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
                        $current_kozak_strength = 3;
                    }
                    elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
                        $current_kozak_strength = 2;
                    }
                }
                if($current_kozak_strength>$uSTOP_KozakStrength){
                    $uSTOP_KozakStrength = $current_kozak_strength;
                    $uSTOP_KozakContext = $current_kozak;
                }

                #if there is an alternative stop codon in the uORF
                if (@stops>1){
                    $uSTOP_AltStop = "True";
                    $uSTOP_AltStopDistance = $stops[1]-$stops[0];
                    $uSTOP_AltStopDistanceToCDS = $length-$stops[1];
                } #if there is no alternative stop codon
                else{
                    $uSTOP_AltStop = "False";
                    $uSTOP_AltStopDistance = "-";
                    $uSTOP_AltStopDistanceToCDS = "-";
                }
            	if (($length-$start_pos) % 3){
                	$uSTOP_FrameWithStart = "outOfFrame";
                    }
                    else{
                    $uSTOP_FrameWithStart = "inFrame";
                    }
            	#find evidence from sorf
            	my $query = "chr".$chr.":".$pos;
            	$uSTOP_evidence=$self->{stop_evidence}->{$query}?"True":"False";
            }
    }

    }

      my %effect = (
    "uSTOP_AltStop" => $uSTOP_AltStop,
    "uSTOP_AltStopDistance" => $uSTOP_AltStopDistance,
    "uSTOP_AltStopDistanceToCDS" => $uSTOP_AltStopDistanceToCDS,
    "uSTOP_FrameWithCDS" => $uSTOP_FrameWithStart,
    "uSTOP_KozakContext" => $uSTOP_KozakContext,
    "uSTOP_KozakStrength" => $kozak_strength{$uSTOP_KozakStrength}? $kozak_strength{$uSTOP_KozakStrength}:$uSTOP_KozakStrength,
    "uSTOP_Evidence" => $uSTOP_evidence,
	);
      my $result=\%effect;

    return $result;
}



    sub creating_ATG{
	#Description: annotate if a five_prime_UTR_variant creates ATG
	#Returntype: hashref

	my ($self, $variant_info,$UTR_info) = @_;
	my %flip;
	$flip{'A'}='T';
	$flip{'C'}='G';
	$flip{'G'}='C';
	$flip{'T'}='A';

	my $chr = $variant_info->{chr};
	my $pos = $variant_info->{pos};
	my $ref = $variant_info->{ref};
	my $alt = $variant_info->{alt};

	my $gene = $UTR_info->{gene};
	my @sorted_starts = @{$UTR_info->{start}};
	my @sorted_ends = @{$UTR_info->{end}};
	my $num_exons = @sorted_starts;
	my @sequence = split //, $UTR_info->{seq};
	my $length = @sequence;
	my $strand = $UTR_info->{strand};

	my @stop_pos = @{$self->get_stopcodon_pos(\@sequence)};

	#return annotators
	my $position_in_ATG = "";
	my $dist_to_start = "";
	my $kozak_strength = "";
	my $conseq = "";
	my $InframeStop = "";
	my $n_ATG_UTR = "";
	my $dist_from_cap = "";
    my $kozak = "";
    my $dist = "";

	#indicate whether the variant creats a ATG
	my $flag;

	#the relative position of input variant in the UTR sequence
	my $n;
	my %chr_position;
	my $utr_position = 0;

	my $result={};

  if ($strand == 1){

  	#create a map from chromosome position to UTR position
  	for (my $m=0; $m<$num_exons; $m++)
  	{
  		for (my $p=$sorted_starts[$m]; $p<=$sorted_ends[$m]; $p++)
  		{
  			$chr_position{$p}=$utr_position;
  			$utr_position++;
  		}
  	}
  	#finding the relative pos of variant in the UTR sequence
  	$n = $chr_position{$pos};

  	#at the 1st position in ATG
  	if(($alt eq 'A') && ($sequence[($n+1)] eq 'T') && ($sequence[($n+2)] eq 'G')){
  		$position_in_ATG=1;
  		if($n>0){
  			$flag=1;
  		}

  	}#at the 2nd position in ATG
    elsif(($n>0) && ($alt eq 'T') && ($sequence[($n-1)] eq 'A') && ($sequence[($n+1)] eq 'G')){
                  $position_in_ATG=2;
                  $flag=1;
    }#at the 3rd position in ATG
  	elsif(($n>1) && ($alt eq 'G') && ($sequence[($n-2)] eq 'A') && ($sequence[($n-1)] eq 'T')){
  		$position_in_ATG=3;
  		$flag=1;
    }

  	if ($flag ==1){

  		#get the pos of A (in ATG) in the UTR sequence
  		my $pos_A = $n-($position_in_ATG-1);

  ################################################################################
  #annotator 1: get the distance to the start codon of the main ORF
  ################################################################################
  		$dist_to_start = $length-$pos_A;
  		$dist_from_cap = $pos_A;
  ################################################################################
  #annotator 2: determine kozak context;
  ################################################################################
  		if ((($pos_A-3)>=0)&&($sequence[($pos_A+3)])){
  			$kozak = $sequence[($pos_A-3)].$sequence[($pos_A-2)].$sequence[$pos_A-1]."ATG".$sequence[$pos_A+3];
  		}
  		else{
  			$kozak = '-';
  		}
  		#get the strength of kozak context
  		if ($kozak !~ /-/){
  			my @split_kozak = split //, $kozak;
  			$kozak_strength = "Weak";
  			if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
  				$kozak_strength = "Strong";
  			}
  			elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
  				$kozak_strength = "Moderate";
  			}
  		}

  ################################################################################
  #annotator 3: Frame with respect the main ORF: Overlapping_inFrame/Overlapping_outofframe
  ################################################################################

  		#if the remainder of the distance to the start codon divided by three is non-zero, then it's Overlapping_outofframe
  		if(($dist_to_start) % 3){
  			$conseq = "outOfFrame";
  		}
  		else{
  			$conseq = "inFrame";
  		}

  ################################################################################
  #annotator 4: Whether there is an inframe stop codon within the 5'UTR
  ################################################################################
		$InframeStop = "False";
  		#get the distance to inframe stop codon
  		
  		foreach my $met (@stop_pos){
  			#if the remainder of the distance to stop codon divided by three is non-zero, then it's out of frame.
  			if (($met-$pos_A) % 3){}
  			#otherwise it's inframe
        #if already finding out a uORF, then finding out the one of min length
  			elsif ((($dist !~/NA/)&&(($met-$pos_A)>0)&&(($met-$pos_A)<$dist))||(($dist =~/NA/)&&(($met-$pos_A)>0)))
  				{
  					$InframeStop = "True";
  					$dist=$met-$pos_A;
  				}

  		}


  		my %effect = (
        "uAUG_KozakContext" => $kozak,
        "uAUG_KozakStrength" => $kozak_strength,
        "uAUG_DistanceToCDS" => $dist_to_start,
        "uAUG_FrameWithCDS" => $conseq,
        "uAUG_InframeStop" => $InframeStop,
        "uAUG_DistanceToInframeStop" => $dist,
        "uAUG_DistanceFrmCap" => $dist_from_cap,
      );
      $result=\%effect;
  	}

  }

  if ($strand == -1){

  	#create a map from chromosome position to UTR position
  	#the exons were arranged in incresing order above

  	for (my $m=$num_exons-1; $m>=0; $m--)
  	{
  		for (my $p=$sorted_ends[$m]; $p>=$sorted_starts[$m]; $p--)
  		{
  			$chr_position{$p}=$utr_position;
  			$utr_position++;
  		}
  	}
  	#finding the relative pos of variant in the sequence (5' to 3' correct strand)
  	$n = $chr_position{$pos};
  	#at the 1st first position in ATG: the reverse nucleotide to A is T.
  	if(($alt eq 'T') && ($sequence[($n+1)] eq 'T') && ($sequence[($n+2)] eq 'G')){
  		$position_in_ATG=1;
  		if($n>0){
  			#the context represented as the forward strand
  			#$trip_context = $flip{$sequence[$n+1]}.$flip{$sequence[$n]}.$flip{$sequence[$n-1]};
  			$flag=1;
  		}
  	}# at the 2nd position in ATG: the complmented nucletide in forward strand is A
  	elsif(($n>0) && ( $alt eq 'A') && ($sequence[($n-1)] eq 'A') && ($sequence[($n+1)] eq 'G')) {
  								$position_in_ATG=2;
                  $flag=1;
    }# at the 3rd position in ATG: the complemented nucleotide in forward strand is C
  	elsif(($n>1) && ($alt eq 'C') && ($sequence[($n-2)] eq 'A') && ($sequence[($n-1)] eq 'T')){
  		$position_in_ATG=3;
  		$flag=1;
    }

  	if ($flag==1){

  		#get the relative pos of A (in ATG) in the UTR sequence
  		my $pos_A = $n-($position_in_ATG-1);

  ################################################################################
  #annotator 1: get the distance to the start codon of the main ORF
  ################################################################################
  		$dist_to_start = $length-$pos_A;
  		$dist_from_cap = $pos_A;
  ################################################################################
  #annotator 2: determine kozak context;
  ################################################################################
  		my $kozak = '';
  		if (($pos_A-3>=0)&&($pos_A+3<=$length-1)) {
  			$kozak = $sequence[($pos_A-3)].$sequence[($pos_A-2)].$sequence[($pos_A-1)]."ATG".$sequence[$pos_A+3];
  		}
  		else{
  			$kozak = '-';
  		}
  		#get the strength of kozak context
  		if ($kozak !~ /-/){
  			my @split_kozak = split //, $kozak;
  			$kozak_strength = "Weak";
  			# for the strong consensus: the -3 and +4 (when the start index of ATG is 1) position match the consensus of kozak sequence;
  			if ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))&&($split_kozak[6] eq 'G')){
  				$kozak_strength = "Strong";
  			}
  			# for moderate consensus: has only 1 of these sites;
  			elsif ((($split_kozak[0] eq 'A')||($split_kozak[0] eq 'G'))||($split_kozak[6] eq 'G')){
  				$kozak_strength = "Moderate";
  			}
  		}

  ################################################################################
  #annotator 3: Frame with respect the main ORF: Overlapping_inFrame/Overlapping_outofframe
  ################################################################################

  		#if the remainder of the distance to the start codon divided by three is non-zero, then it's Overlapping_outofframe
  		if(($dist_to_start) % 3){
  			$conseq = "outOfFrame";
  		}
  		else{
  			$conseq = "inFrame";
  		}

  ################################################################################
  #annotator 4: Whether there is an inframe stop codon within the 5'UTR
  ################################################################################

		$InframeStop = "False";
  		#get the distance to inframe stop codon
  		my $dist = 'NA';
  		foreach my $met (@stop_pos){
  			#if the remainder of the distance to stop codon divided by three is non-zero, then it's out of frame.
  			if (($met-$pos_A) % 3){}
  			#otherwise it's inframe
  			elsif ((($dist !~/NA/)&&(($met-$pos_A)>0)&&(($met-$pos_A)<$dist))||(($dist =~/NA/)&&(($met-$pos_A)>0)))
  				{
  					$InframeStop = "True";
  					$dist=$met-$pos_A;
  				}

  		}

  my %effect = (
    "uAUG_KozakContext" => $kozak,
    "uAUG_KozakStrength" => $kozak_strength,
    "uAUG_DistanceToCDS" => $dist_to_start,
    "uAUG_FrameWithCDS" => $conseq,
    "uAUG_InframeStop" => $InframeStop,
    "uAUG_DistanceToInframeStop" => $dist,
    "uAUG_DistanceFrmCap" => $dist_from_cap,
  );
      $result=\%effect;
  	}

  }
   return $result;

}

    sub count_number_ATG{

      #Description: count the number of existing ATGs in the five prime UTR sequence
      #Return: two numbers: the number of uORFs and the number of oORFs
      #Returntype: listref

      my ($self,$seq) = @_;
      my @sequence = split //, $seq;
	  my $length = @sequence;

      my @atg_pos = @{$self->get_ATG_pos(\@sequence)};
      my @mes_pos = @{$self->get_stopcodon_pos(\@sequence)};

      my $atg_num = $self->get_ATG_pos(\@sequence);
      my $inframe_stop_num=0;
      my $outofframe_atg_num=0;
      my $inframeORF_num=0;
      foreach my $atg (@atg_pos){
        my $flag=0; #indicate whether there is a stop codon with respect to this ATG
        foreach my $mes (@mes_pos){
		        if ((($mes-$atg) % 3)||($mes-$atg)<0){}
		        else{
            $flag=1;
            last;
          }
        }
        #if there is no stop codon, then look at whether it's Out_of_frame or Inframe
        if($flag==1){$inframe_stop_num++;}
		else{
         (($length-$atg) % 3 !=0)?$outofframe_atg_num++:$inframeORF_num++;
       }
      }

      my %existing_uORF = (
      "existing_uORFs" => $inframe_stop_num,
      "existing_oORFs" => $outofframe_atg_num,
      "existing_inframeORFs" => $inframeORF_num,);

        return \%existing_uORF;

    }

    sub existing_uORF{

        #Description: obtaining the relative coordinates of start and end pos of existing uORF in the five prime UTR sequence
        #Return: hashref(key:the pos of the first nucleotidie of start codon; value:all the positions of the first nucleotide of stop codon)

        my ($self,$seq) = @_;
        my $length = @{$seq};

        my @atg_pos = @{$self->get_ATG_pos($seq)};
        my @mes_pos = @{$self->get_stopcodon_pos($seq)};
        my %uORF;

        foreach my $atg (@atg_pos){
            my @end_pos;
            foreach my $mes (@mes_pos){
                if ((($mes-$atg) % 3)||($mes-$atg)<0){}
                else{
                    push @end_pos,$mes;
                }
            }
            if (@end_pos){
            $uORF{$atg}=[@end_pos];
        }

    }
        return \%uORF;
    }


    sub get_ATG_pos{

      #Description: get all the relative position of ATG in the five prime UTR sequence
      #Returntype: listref

      my ($self,$seq) = @_;

      my @sequence = @{$seq};
      my $length = @sequence;

      my @atg_pos;
      for (my $seq_n=0;$seq_n<=$length-3;$seq_n++){
        if(($sequence[($seq_n)] eq 'A')&&($sequence[($seq_n+1)] eq 'T')&&($sequence[($seq_n+2)] eq 'G')){
          push @atg_pos,$seq_n;
        }
      }
      return \@atg_pos;

    }

    sub get_stopcodon_pos{
	#Description: get all the relative position of stop codons in the five prime UTR sequence
	#Returntype: listref
	my ($self,$seq) = @_;

  	my @sequence = @{$seq};
  	my $length = @sequence;

	my @met_pos;
	for (my $seq_n=0; $seq_n<$length; $seq_n++){
	if ((($sequence[$seq_n] eq 'T')&&($sequence[$seq_n+1] eq 'A')&&($sequence[$seq_n+2] eq 'A'))
	||(($sequence[$seq_n] eq 'T')&&($sequence[$seq_n+1] eq 'A')&&($sequence[$seq_n+2] eq 'G'))
	||(($sequence[$seq_n] eq 'T')&&($sequence[$seq_n+1] eq 'G')&&($sequence[$seq_n+2] eq 'A')))
				{
					push @met_pos,$seq_n;
				}
	}
	#return a reference
	return \@met_pos;
}


 1;
