#!/usr/bin/perl -w
use strict;

my $min_length = 19;
# Specifications
my $line_counter = 0; my $f_counter = 0; my $gatc_counter; my %hash; my $flag; my $k;
my @line_labels = qw(name seq spare qscore);

my $filename = $ARGV[0];
unless ($filename =~ /(.*)\.(fastq|fq)/) {die"filename does not match fq"};
my $filename_out= $1."_REdig.$2"; my $filename_dump = $1."_merge_dump.$2";
my $file_name=$1;
my $file_path="";
if ($file_name =~ /(.*\/)(\V++)/) {$file_path = $1; $file_name = $2};

unless (open(FH, $filename)) {print "Cannot open file $filename\n"; exit;}

# opens a file in append modes for the output of the data
open FHOUT, ">$filename_out" or die $!;   
  

while ($hash{$line_labels[$f_counter]}=<FH>)  #assigns each fq line to the hash in batches of 4
{
chomp $hash{$line_labels[$f_counter]};
$f_counter++; $line_counter++;


if ($f_counter==4)
    {
    # name formats @HISEQ2000:376:C2399ACXX:8:1101:1749:1893 1:N:0:GAGTTAGT run1
    # name formats @HISEQ2000:376:C2399ACXX:8:1101:1749:1893 2:N:0:GAGTTAGT run2
    # /(.*):(.*):(.*):(.*):(.*):(.*):(\d++) (\d):(.*):(.*):(.*)/
    if ($hash{"name"} =~ /(.*:.*:.*:.*:.*:.*:\d++) (\d):(.*:.*:.*)/)
        {
        $hash{"PE"}=$2;
        $hash{"new name"} = $1;
        #$hash{"new name end"} = " ".$2.":".$3
        }
    
    if ($hash{"seq"} =~ /GATC/)
        {

         my @gatc_splits = split/GATC/, $hash{"seq"};

             for (my $i=0; $i<$#gatc_splits+1;$i++)
             {
             if ($i==0) #first fragment
                {
                $hash{"split$i"}{"sseq"}= "$gatc_splits[$i]GATC";
                $hash{"split$i"}{"sqscore"}= substr ($hash{qscore},0,length $hash{"split$i"}{"sseq"});
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i;
                #$hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.$hash{"new name end"};
                }
             if ($i!=0 and $i != $#gatc_splits) #middle fragment
                {
                $hash{"split$i"}{"sseq"}= "GATC$gatc_splits[$i]GATC";
                my $offset=0;
                for (my$j=0; $j<$i; $j++){$offset = $offset -4 + length $hash{"split$j"}{"sseq"};} # calculates the offset by looping through the lengths of the left hand side fragments and summing them
                $hash{"split$i"}{"sqscore"}= substr ($hash{qscore},$offset,length $hash{"split$i"}{"sseq"});
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i;
                #$hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.$hash{"new name end"};                
                }
             if ($i==$#gatc_splits and $i!=0) #last fragment if there is more than one fragment
                {
                $hash{"split$i"}{"sseq"}= "GATC$gatc_splits[$i]";
                $hash{"split$i"}{"sqscore"}= substr ($hash{qscore},-length $hash{"split$i"}{"sseq"},length $hash{"split$i"}{"sseq"});
                $hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i;
                #$hash{"split$i"}{"sname"}= $hash{"new name"}.":PE".$hash{"PE"}.":".$i.$hash{"new name end"};                   
                }         
         
             if (length $hash{"split$i"}{"sseq"}>$min_length)
                {
                print FHOUT $hash{"split$i"}{"sname"}."\n".$hash{"split$i"}{"sseq"}."\n+\n".$hash{"split$i"}{"sqscore"}."\n";   
                }
             }

         $gatc_counter++;
        }
    else
    {
             print FHOUT $hash{"new name"}.":PE".$hash{"PE"}.":0"."\n";
             #print FHOUT $hash{"new name"}.":PE".$hash{"PE"}.":0".$hash{"new name end"}."\n";
             print FHOUT $hash{"seq"}."\n";  #error check ."\t".length $hash{"split$i"}{"sseq"};
             print FHOUT "+\n";
             print FHOUT $hash{"qscore"}."\n";   #error check "\t".length $hash{"split$i"}{"sqscore"};        
    }
        
    #prints the data in the hash: for (my $i=0; $i<4; $i++){print $i.$hash{$line_labels[$i]}."\n"}print "\n";

    $f_counter=0

    }

#if ($line_counter>1000000){print "$line_counter lines reviewed\n$gatc_counter DpnII sites found\n";exit #}   
}
print "dpnII.pl command run on file: $filename\n$line_counter lines reviewed\n$gatc_counter DpnII sites found\n"