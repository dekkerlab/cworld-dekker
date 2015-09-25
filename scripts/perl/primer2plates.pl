#!/usr/bin/perl -w

use English;
use warnings;
use strict;
use Carp qw(carp cluck croak confess);
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use POSIX qw(ceil floor);
use List::Util qw[min max];
use Cwd 'abs_path';
use Cwd;

use cworld::dekker;

my $tool=(split(/\//,abs_path($0)))[-1];

sub check_options {
    my $opts = shift;

    my ($inputFile,$verbose,$output);
    
    my $ret={};
    
    if( $opts->{ inputFile } ) {
        $inputFile = $opts->{ inputFile };
    } else {
        print STDERR "\nERROR: Option inputFile|i is required.\n";
        help();
    }
    
    if( exists($opts->{ verbose }) ) {
        $verbose = $opts->{ verbose };
    } else {
        $verbose = 0;
    }
    
    if( exists($opts->{ output }) ) {
        $output = $opts->{ output };
    } else {
        $output = "myPrimers";
    }
   
   $ret->{ inputFile }=$inputFile;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    
    return($ret,$inputFile,$verbose,$output);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\tlayout primers in 96-well plate format\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl primer2plates.pl [OPTIONS] -i <inputPrimers>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file, multiple files allowed");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");
    
    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script can layout primers in 96-well plate format
    Matrix should have row/col headers in the standard my5C format.
    Matrix can be TXT or gzipped TXT.
    See website for matrix format details.\n";
    
    print STDERR "\n";
    
    print STDERR "Contact:
    Bryan R. Lajoie
    Dekker Lab 2015
    https://github.com/blajoie/cworld-dekker
    http://my5C.umassmed.edu";
    
    print STDERR "\n";
    print STDERR "\n";
    
    exit;
}

my %options;
my $results = GetOptions( \%options,'inputFile|i=s','verbose|v','output|o=s') or croak help();
my ($ret,$inputFile,$verbose,$output) = &check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);

my %well2row=();

$well2row{0}="A";
$well2row{1}="B";
$well2row{2}="C";
$well2row{3}="D";
$well2row{4}="E";
$well2row{5}="F";
$well2row{6}="G";
$well2row{7}="H";

my %wellInc=();
$wellInc{ FOR }=0;
$wellInc{ REV }=0;
$wellInc{ LFOR }=0;
$wellInc{ LREV }=0;
my %plateInc=();
$plateInc{ FOR }=0;
$plateInc{ REV }=0;
$plateInc{ LFOR }=0;
$plateInc{ LREV }=0;

open(OUTF,outputWrapper($output."_F.txt")) or croak "Could not open file [$output] - $!";
open(OUTR,outputWrapper($output."_R.txt")) or croak "Could not open file [$output] - $!";
open(OUTLF,outputWrapper($output."_LF.txt")) or croak "Could not open file [$output] - $!";
open(OUTLR,outputWrapper($output."_LR.txt")) or croak "Could not open file [$output] - $!";

my $lastPrimerSetName="";
my $lastPrimerType="";

my $found_flag=0;
my $lineNum=0;
my %col2header=();

print STDERR "converting primers to plate ...\n" if($verbose);

open(IN,inputWrapper($inputFile)) or croak "Could not open file [$inputFile] - $!";
while(my $line = <IN>) {
    chomp($line);
    next if(($line eq "") or ($line =~ m/^#/));

    if($lineNum == 0) {
        
        my @tmp=split(/\t/,$line);
        for(my $i=0;$i<@tmp;$i++) {
            $col2header{$tmp[$i]}=$i;
        }

        croak "could not find the 'PRIMER_NAME' column...\n" if(!exists($col2header{ PRIMER_NAME }));
        croak "could not find the 'PRIMER_SEQUENCE' column...\n" if(!exists($col2header{ PRIMER_SEQUENCE }));
        
    } else {
        
        my @tmp=split(/\t/,$line);
        my $primerName=$tmp[$col2header{ PRIMER_NAME }];
        my @pNameArr=split(/\|/,$primerName);
        $primerName=$pNameArr[0];
        
        my $primerSequence=$tmp[$col2header{ PRIMER_SEQUENCE }];
        
        my @name_tmp=split(/_/,$primerName);
        my $primerSetName=$name_tmp[2].".".$name_tmp[3];
        my $primerType=$name_tmp[3];
        
        if($lastPrimerSetName eq "") { $lastPrimerSetName = $primerSetName; }
        
        if($primerSetName ne $lastPrimerSetName) {
            
            if($found_flag) {
                #add 3 gaps in the plate
                for(my $i=0;$i<3;$i++) {
                    
                    my $wellNum=$wellInc{$lastPrimerType};
                    my $plateNum=$plateInc{$lastPrimerType};
                    
                    my $row=$well2row{(floor($wellNum/12))};
                    my $col=($wellNum%12)+1;
                    
                    my $plateName="";
                    if($lastPrimerType eq "FOR") {
                        $plateName=$output."_F_".$plateNum;
                        print OUTF "$plateName\t$row\t$col\tEMPTY\t\tEMPTY\n";
                        $found_flag=1;
                    } elsif($lastPrimerType eq "REV") {
                        $plateName=$output."_R_".$plateNum;
                        print OUTR "$plateName\t$row\t$col\tEMPTY\t\tEMPTY\n";
                        $found_flag=1;
                    } elsif($lastPrimerType eq "LFOR") {
                        $plateName=$output."_LF_".$plateNum;
                        print OUTLF "$plateName\t$row\t$col\tEMPTY\t\tEMPTY\n";
                        $found_flag=1;
                    } elsif($lastPrimerType eq "LREV") {
                        $plateName=$output."_LR_".$plateNum;
                        print OUTLR "$plateName\t$row\t$col\tEMPTY\t\tEMPTY\n";
                        $found_flag=1;
                    } else {
                        print STDERR "found an invalid primerType($primerType)...\n";
                        exit;
                    }
                    
                    $wellInc{$lastPrimerType}++;
                    #if > 95 - make a new plate and reset wellInc
                    if($wellInc{$lastPrimerType} > 95) { 
                        $wellInc{$lastPrimerType}=0; 
                        $plateInc{$lastPrimerType}++;
                    }
                    
                }
                $found_flag=0;
            }
        }
        $lastPrimerSetName=$primerSetName;
        $lastPrimerType=$primerType;
        
        my $wellNum=$wellInc{$primerType};
        my $plateNum=$plateInc{$primerType};
        
        my $row=$well2row{(floor($wellNum/12))};
        my $col=($wellNum%12)+1;
            
        my $plateName="";
        if($primerType eq "FOR") {
            $plateName=$output."_F_".$plateNum;
            print OUTF "$plateName\t$row\t$col\t$primerName\t\t$primerSequence\n";
            $found_flag=1;
        } elsif($primerType eq "REV") {
            $plateName=$output."_R_".$plateNum;
            print OUTR "$plateName\t$row\t$col\t$primerName\t\t$primerSequence\n";
            $found_flag=1;
        } elsif($primerType eq "LFOR") {
            $plateName=$output."_LF_".$plateNum;
            print OUTLF "$plateName\t$row\t$col\t$primerName\t\t$primerSequence\n";
            $found_flag=1;
        } elsif($primerType eq "LREV") {
            $plateName=$output."_LR_".$plateNum;
            print OUTLR "$plateName\t$row\t$col\t$primerName\t\t$primerSequence\n";
            $found_flag=1;
        } else {
            print STDERR "found an invalid primerType($primerType)...\n";
            exit;
        }
        
        $wellInc{$primerType}++;
        #if > 95 - make a new plate and reset wellInc
        if($wellInc{$primerType} > 95) { 
            $wellInc{$primerType}=0; 
            $plateInc{$primerType}++;
        }
    }
    $lineNum++;
}

print STDERR "\tdone\n" if($verbose);

close(OUTF);
close(OUTR);

print STDERR "\n" if($verbose);
