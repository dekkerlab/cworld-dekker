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

    my ($inputMatrix,$verbose,$output);
    
    my $ret={};
    
    if( exists($opts->{ inputMatrix }) ) {
        $inputMatrix = $opts->{ inputMatrix };
    } else {
        print STDERR "\nERROR: Option inputMatrix|i is required.\n";
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
        $output = "";
    }
    
    $ret->{ inputMatrix }=$inputMatrix;
    $ret->{ verbose }=$verbose;
    $ret->{ output }=$output;
    
    return($ret,$inputMatrix,$verbose,$output);
}

sub intro() {
    print STDERR "\n";
    
    print STDERR "Tool:\t\t".$tool."\n";
    print STDERR "Version:\t".$cworld::dekker::VERSION."\n";
    print STDERR "Summary:\ttransform symmetrical matrix into axis-seperated [5C headers only]\n";
    
    print STDERR "\n";
}

sub help() {
    intro();
    
    print STDERR "Usage: perl symmetrical2seperate.pl [OPTIONS] -i <inputMatrix>\n";
    
    print STDERR "\n";
    
    print STDERR "Required:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-i", "[]", "input matrix file");
    
    print STDERR "\n";
    
    print STDERR "Options:\n";
    printf STDERR ("\t%-10s %-10s %-10s\n", "-v", "[]", "FLAG, verbose mode");
    printf STDERR ("\t%-10s %-10s %-10s\n", "-o", "[]", "prefix for output file(s)");

    print STDERR "\n";
    
    print STDERR "Notes:";
    print STDERR "
    This script turns any symmetrical matrix into a axis-seperated matrix.  FORWARD on y-axis, REVERSE on x-axis.
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

sub seperateAxes($$) {
    my $inc2header=shift;
    my $header2inc=shift;
    
    my $numTotalHeaders=keys(%{$header2inc->{ xy }});
    
    my $seperated_inc2header={};
    my $seperated_header2inc={};
    
    my $for=0;
    my $rev=0;
    
    for(my $xy=0;$xy<$numTotalHeaders;$xy++) {
        my $header=$inc2header->{ xy }->{$xy};
        my $headerObject=getHeaderObject($header);
        my $headerPrimerType=$headerObject->{ primerType };
        if(($headerPrimerType eq "FOR") or ($headerPrimerType eq "LFOR")) {
            $seperated_inc2header->{ y }->{$for}=$header;
            $seperated_header2inc->{ y }->{$header}=$for;
            $for++;
        }
        if(($headerPrimerType eq "REV") or ($headerPrimerType eq "LREV")) {
            $seperated_inc2header->{ x }->{$rev}=$header;
            $seperated_header2inc->{ x }->{$header}=$rev;
            $rev++;
        }    
    }
    
    croak "invalid headers - can only seperate my5C formatted headers (5C_2410_EMS03_REV_898|hg18|chr7:118385397-118394195)" if($rev == 0);
    croak "invalid headers - can only seperate my5C formatted headers (5C_2410_EMS03_FOR_900|hg18|chr7:118399267-118401405)" if($for == 0);
    
    return($seperated_inc2header,$seperated_header2inc);
}

sub desymmetricizeData($$$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $seperated_inc2header=shift;
    my $seperated_header2inc=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $missingValue=$matrixObject->{ missingValue };
    
    my $numYHeaders=keys(%{$seperated_header2inc->{ y }});
    my $numXHeaders=keys(%{$seperated_header2inc->{ x }});
    
    my $seperatedMatrix={};
    
    for(my $y=0;$y<$numYHeaders;$y++) {
    
        my $yHeader=$seperated_inc2header->{ y }->{$y};
        my $symmetrical_y_index=-1;
        $symmetrical_y_index=$header2inc->{ y }->{$yHeader};
    
        for(my $x=0;$x<$numXHeaders;$x++) {
            
            my $xHeader=$seperated_inc2header->{ x }->{$x};
            my $symmetrical_x_index=-1;
            $symmetrical_x_index=$header2inc->{ x }->{$xHeader};
        
            my $cScore=$matrixObject->{ missingValue };
            
            $cScore=$matrix->{$symmetrical_y_index}->{$symmetrical_x_index} if(defined($matrix->{$symmetrical_y_index}->{$symmetrical_x_index}));
            $seperatedMatrix->{$y}->{$x}=$cScore;
            
        }
    }
    
    return($seperatedMatrix);
}

my %options;
my $results = GetOptions( \%options,'inputMatrix|i=s','verbose|v','output|o=s') or croak help();
my ($ret,$inputMatrix,$verbose,$output)=check_options( \%options );

intro() if($verbose);

#get the absolute path of the script being used
my $cwd = getcwd();
my $fullScriptPath=abs_path($0);
my @fullScriptPathArr=split(/\//,$fullScriptPath);
@fullScriptPathArr=@fullScriptPathArr[0..@fullScriptPathArr-3];
my $scriptPath=join("/",@fullScriptPathArr);
my $commentLine=getScriptOpts($ret,$tool);

croak "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

# get matrix information
my $matrixObject=getMatrixObject($inputMatrix,$output,$verbose);
my $inc2header=$matrixObject->{ inc2header };
my $header2inc=$matrixObject->{ header2inc };
my $numYHeaders=$matrixObject->{ numYHeaders };
my $numXHeaders=$matrixObject->{ numXHeaders };
my $missingValue=$matrixObject->{ missingValue };
my $symmetrical=$matrixObject->{ symmetrical };
my $equalHeaderFlag=$matrixObject->{ equalHeaderFlag };
my $inputMatrixName=$matrixObject->{ inputMatrixName };
$output=$matrixObject->{ output };

my $matrix={};
($matrix)=getData($inputMatrix,$matrixObject,$verbose);
print STDERR "\tdone\n" if($verbose);

print STDERR "sepeating axes ...\n" if($verbose);
my ($seperated_inc2header,$seperated_header2inc)=seperateAxes($inc2header,$header2inc);
$numYHeaders=keys(%{$seperated_header2inc->{ y }});
$numXHeaders=keys(%{$seperated_header2inc->{ x }});
print STDERR "\tnumYHeaders\t$numYHeaders\n" if($verbose);
print STDERR "\tnumxHeaders\t$numXHeaders\n" if($verbose);

print STDERR "\n" if($verbose);

print STDERR "de-symmetriciz'n data...\n" if($verbose);
$matrix=desymmetricizeData($matrixObject,$matrix,$seperated_inc2header,$seperated_header2inc);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);

my $symmetricalMatrixFile = $output.".nonsymmetrical.matrix.gz";
print STDERR "writing matrix (symmetricalMatrixFile) ...\n" if($verbose);
writeMatrix($matrix,$seperated_inc2header,$symmetricalMatrixFile,$matrixObject->{ missingValue },$commentLine);
print STDERR "\tdone\n" if($verbose);

print STDERR "\n" if($verbose);