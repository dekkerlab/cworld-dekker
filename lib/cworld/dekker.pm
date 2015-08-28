#!/usr/bin/perl -w

package cworld::dekker;

use 5.006;
use strict;
use warnings;
use Carp qw(carp cluck confess confess);
use POSIX qw(ceil floor strftime);
use List::Util qw[min max];
use Cwd 'abs_path';
use Cwd;

=head1 NAME

cworld::dekker - perl module and collection of utility/analysis scripts for C data (3C, 4C, 5C, Hi-C)

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Large collection of functions and utility scripts for processing of C data

    use cworld::dekker;

=head1 AUTHOR

Bryan R. Lajoie, C<< <bryan.lajoie at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-cworld-dekker at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=cworld-dekker>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc cworld::dekker

You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=cworld-dekker>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/cworld-dekker>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/cworld-dekker>

=item * Search CPAN

L<http://search.cpan.org/dist/cworld-dekker/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2015 Bryan R. Lajoie.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=cut

require Exporter;

our @ISA = qw(Exporter);

our @EXPORT = qw(autoScale autoSize badFormat baseName
             calculateColorPalette calculateLoess calculateLog2Ratio
             calculateObsMinusExp calculateTransExpected calculateZscore
             classifyInteraction classifyInteractionDistance
             combineBedFiles commify compareMatrices
             correctMatrix correlateMatrices createTmpDir deGroupHeader
             flipBool getAvailableColors getColorIndex
             getComputeResource getData getDate getFileName
             getFilePath getHeaderSpacing getMatrixAttributes getMatrixObject
             getMatrixSum getNARows getNumberOfLines getRowColFactor
             getHeaderObject getRestrictionEnzymeSequences
             getRowSum getShortFileName getSmallUniqueString
             getUniqueString getUserHomeDirectory header2subMatrix headers2bed
             initHeatmap inputWrapper intersectBED intersectHeaders isOverlapping
             isSymmetrical listStats loadBED logTransformMatrix matrix2distance
             matrix2inputlist matrix2listfile matrix2pairwise midpointBedFile
             normalizeMatrix outputWrapper parseHeaders processMatrixFile
             readLoessFile removeDiagonal removeFileExtension removeTmpDir round
             roundNearest splitCoordinate stitchMatrices stripChromosomeGroup 
             translateFlag transposeMatrix getInteractionDistance updateMatrixObject validateBED
             validateIdenticalMatrixStructure validateLoessObject validateZoomCoordinate
             writeMatrix);
             
our @EXPORT_OK = @EXPORT;

=head2 getDate

 Title     : getDate
 Usage     : $time=getDate()
 Function  : returns current date/time.
 Returns   : string
 Argument  : None

=cut

sub getDate() {
    my $time = strftime '%I:%M:%S %P, %m/%d/%Y', localtime;
    
    return($time);
}

=head2 commify

 Title     : commify
 Usage     : $num_string=commify($num)
 Function  : adds commas to numbers (every 1000)
 Returns   : string
 Argument  : number to format

=cut

sub commify {
   (my $num = shift) =~ s/\G(\d{1,3})(?=(?:\d\d\d)+(?:\.|$))/$1,/g; 
   return $num; 
}
    
=head2 writeMatrix

 Title     : writeMatrix
 Usage     : writeMatrix(...)
 Function  : write a text (tsv) matrix file
 Returns   : None
 Argument  : 2D hash, row/col hash, file path

=cut

sub writeMatrix($$$;$$) {
    #required
    my $matrix=shift;
    my $inc2header=shift;
    my $matrixFile=shift;
    #optional
    my $missingValue=0;
    $missingValue=shift if @_;
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    my $numYHeaders=keys(%{$inc2header->{ y }});
    my $numXHeaders=keys(%{$inc2header->{ x }});
    
    open(OUT,outputWrapper($matrixFile)) or confess "Could not open file [$matrixFile] - $!";
    
    print OUT $numXHeaders."x".$numYHeaders;
    
    for(my $x=0;$x<$numXHeaders;$x++) {
        my $xHeader=$inc2header->{ x }->{$x};
        print OUT "\t".$xHeader;
    }
    
    print OUT "\n";
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        
        print OUT "$yHeader";
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            
            my $cScore=$missingValue;
            $cScore = $matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            $cScore = "NA" if(($cScore eq "") or ($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i) or ($cScore eq "NA"));
            $cScore = sprintf "%.".$sigDigits."f", $cScore if(($cScore ne "NA") and ($cScore !~ /^[+-]?\d+$/));
            
            print OUT "\t$cScore";
            
        }
        
        print OUT "\n" if($y != ($numYHeaders-1))
    }
    
    close(OUT);

}

=head2 calculateZscore

 Title     : calculateZscore
 Usage     : $zScoreMatrix=calculateZscore(...)
 Function  : translate matrix into z-score matrix
 Returns   : 2D hash 
 Argument  : matrixObject, 2D hash, loess hash

=cut

sub calculateZscore($$$;$$$) {
    #required
    my $matrixObject=shift;
    my $matrix=shift;
    my $loess=shift;
    #optional
    my $cisApproximateFactor=1;
    $cisApproximateFactor=shift if @_;
    my $excludeZero=0;
    $excludeZero=shift if @_;
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $verbose=$matrixObject->{ verbose };
    my $output=$matrixObject->{ output };

    my $zScoreMatrix={};
    
    my ($y,$x);
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);
        
            my $cScore=$matrixObject->{ missingValue };
            $cScore=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
                
            $cScore = "NA" if($cScore eq ""); 
            $cScore = "NA" if(($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i));
            my $zScore = "NA";
            
            if( ($cScore ne "NA") and (exists($loess->{$interactionDistance}->{ loess })) ) {
            
                confess "distance ($interactionDistance) does not exist!" if(!exists($loess->{$interactionDistance}));
                
                my $loessValue=$loess->{$interactionDistance}->{ loess };
                my $loessStdev=$loess->{$interactionDistance}->{ stdev };
                
                $zScore = (($cScore-$loessValue)/$loessStdev) if(($loessValue ne "NA") and ($loessStdev != 0));
            }
            $zScore="NA" if(($excludeZero) and ($cScore == 0));
            
            # truncate numbers to minimal digits
            $zScore = sprintf "%.".$sigDigits."f", $zScore if($zScore ne "NA");
            
            $zScoreMatrix->{$y}->{$x}=$zScore if($zScore ne "NA");
        }
        my $pcComplete = round((($y/($numYHeaders-1))*100),2);
        print STDERR "\e[A" if(($verbose) and ($y != 0));
        printf STDERR "\t%.2f%% complete (".$y."/".($numYHeaders-1).")...\n", $pcComplete if($verbose);
    }
    
    return($zScoreMatrix);

}

=head2 calculateLog2Ratio

 Title     : calculateLog2Ratio
 Usage     : $log2ratioMatrix=calculateLog2Ratio(...)
 Function  : translate matrix into log2(obs/exp) matrix
 Returns   : 2D hash 
 Argument  : matrixObject, 2D hash, loess hash

=cut

sub calculateLog2Ratio($$$;$$$) {
    #required
    my $matrixObject=shift;
    my $matrix=shift;
    my $loess=shift;
    #optional
    my $cisApproximateFactor=1;
    $cisApproximateFactor=shift if @_;
    my $excludeZero=0;
    $excludeZero=shift if @_;
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $verbose=$matrixObject->{ verbose };
    
    my $log2ratioMatrix={};
    
    my ($y,$x);
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);
        
            my $cScore=$matrixObject->{ missingValue };
            $cScore=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
                        
            $cScore = "NA" if($cScore eq ""); 
            $cScore = "NA" if(($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i));
            my $log2ratio = "NA";
            
            if( ($cScore ne "NA") and (exists($loess->{$interactionDistance}->{ loess })) ){
            
                confess "distance ($interactionDistance) does not exist!" if(!exists($loess->{$interactionDistance}));
                
                my $loessValue=$loess->{$interactionDistance}->{ loess };
                my $loessStdev=$loess->{$interactionDistance}->{ stdev };
                
                $log2ratio = log($cScore/$loessValue)/log(2) if(($loessValue ne "NA") and ($loessValue > 0) and ($cScore > 0));
            }
            $log2ratio="NA" if(($excludeZero) and ($cScore == 0));
            
            # truncate numbers to minimal digits
            $log2ratio = sprintf "%.".$sigDigits."f", $log2ratio if($log2ratio ne "NA");
            
            $log2ratioMatrix->{$y}->{$x}=$log2ratio if($log2ratio ne "NA");
        }
        my $pcComplete = round((($y/($numYHeaders-1))*100),2);
        print STDERR "\e[A" if(($verbose) and ($y != 0));
        printf STDERR "\t%.2f%% complete (".$y."/".($numYHeaders-1).")...\n", $pcComplete if($verbose);
    }
    
    return($log2ratioMatrix);

}

=head2 calculateObsMinusExp

 Title     : calculateObsMinusExp
 Usage     : $log2ratioMatrix=calculateObsMinusExp(...)
 Function  : translate matrix into obs-exp matrix
 Returns   : 2D hash 
 Argument  : matrixObject, 2D hash, loess hash

=cut

sub calculateObsMinusExp($$$;$$$) {
    #required
    my $matrixObject=shift;
    my $matrix=shift;
    my $loess=shift;
    #optional
    my $cisApproximateFactor=1;
    $cisApproximateFactor=shift if @_;
    my $excludeZero=0;
    $excludeZero=shift if @_;
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $verbose=$matrixObject->{ verbose };
    
    my $obsMinusExpMatrix={};
    
    my ($y,$x);
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);
        
            my $cScore=$matrixObject->{ missingValue };
            $cScore=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
                        
            $cScore = "NA" if($cScore eq ""); 
            $cScore = "NA" if(($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i));
            
            my $obsMinusExp = "NA";
            if( ($cScore ne "NA") and (exists($loess->{$interactionDistance}->{ loess })) ){
            
                confess "distance ($interactionDistance) does not exist!" if(!exists($loess->{$interactionDistance}));
                
                my $loessValue=$loess->{$interactionDistance}->{ loess };
                my $loessStdev=$loess->{$interactionDistance}->{ stdev };
                
                $obsMinusExp = $cScore-$loessValue if($loessValue ne "NA");
            }
            $obsMinusExp="NA" if(($excludeZero) and ($cScore == 0));
                
            # truncate numbers to minimal digits
            $obsMinusExp = sprintf "%.".$sigDigits."f", $obsMinusExp if($obsMinusExp ne "NA");
            
            $obsMinusExpMatrix->{$y}->{$x}=$obsMinusExp if($obsMinusExp ne "NA");
        }
        my $pcComplete = round((($y/($numYHeaders-1))*100),2);
        print STDERR "\e[A" if(($verbose) and ($y != 0));
        printf STDERR "\t%.2f%% complete (".$y."/".($numYHeaders-1).")...\n", $pcComplete if($verbose);
    }
    
    return($obsMinusExpMatrix);

}

=head2 listStats

 Title     : listStats
 Usage     : $dataHash=listStats(...)
 Function  : return math stats on input arr
 Returns   : hash 
 Argument  : array ref, trim pc[0-1]

=cut

sub listStats($;$) {
    # required
    my $listRef=shift;
    # optional
    my $trimPC=0;
    $trimPC=shift if @_;
    
    #flip trim for top/bottom
    $trimPC = ($trimPC / 2);
    
    confess "invalid trimPC value (0-1 range)" if(($trimPC > 1) or ($trimPC < 0));
    
    my $listArrSize=0;
    $listArrSize=@$listRef if(defined($listRef));
    
    my ($mean,$stdev,$variance,$median,$q1,$q3,$iqr,$min,$max,$count,$sum,$binary,$geomean);
    $mean=$stdev=$variance=$median=$q1=$q3=$iqr=$min=$max=$count=$sum=$binary=$geomean=0;
    
    my $skipGeoMean=0;
    my $iqrMean="NA";
    my $zeroPC=0;
    
    # listArr and sort list
    my @listArr=@$listRef;
    @listArr = sort { $a <=> $b } @listArr;
    
    # validate non null list
    for(my $i=0;$i<$listArrSize;$i++) {
        my $tmpVal=$listArr[$i];
        confess "found null in listStats(list)" if($tmpVal eq "NA");
    }
    
    confess "the list is EMPTY!" if($listArrSize == 0);
    
    #quartiles
    my $medianPosition=0.5;
    my $q1Position=0.25;
    my $q3Position=0.75;
    
    if(($listArrSize % 2)) {
        #median
        my $medianIndex = floor($listArrSize * $medianPosition);
        $median = $listArr[$medianIndex];
        #quartile 1
        my $q1Index = floor($listArrSize * $q1Position);
        $q1 = $listArr[$q1Index];
        #quartile 3
        my $q3Index = floor($listArrSize * $q3Position);
        $q3 = $listArr[$q3Index];
        
    } else {
        #median
        my $medianIndexLeft = floor($listArrSize * $medianPosition) - 1;
        my $medianIndexRight = $medianIndexLeft + 1;        
        my $medianLeft = $listArr[$medianIndexLeft];
        my $medianRight = $listArr[$medianIndexRight];
        $median = (($medianLeft + $medianRight) / 2);        
        #quartile 1
        my $q1IndexLeft = floor($listArrSize * $q1Position) - 1;
        my $q1IndexRight = $q1IndexLeft + 1;        
        my $q1Left = $listArr[$q1IndexLeft];
        my $q1Right = $listArr[$q1IndexRight];
        $q1 = (($q1Left + $q1Right) / 2);        
        #quartile 3
        my $q3IndexLeft = floor($listArrSize * $q3Position) - 1;
        my $q3IndexRight = $q3IndexLeft + 1;        
        my $q3Left = $listArr[$q3IndexLeft];
        my $q3Right = $listArr[$q3IndexRight];
        $q3 = (($q3Left + $q3Right) / 2);        
    }
    
    $iqr = ($q3-$q1);
    my $iqrLowerBound = $q1-(1.5*$iqr);
    my $iqrUpperBound = $q3+(1.5*$iqr);
    my $iqrSum=0;
    my $iqrCount=0;
    for(my $i=0;$i<$listArrSize;$i++) {
        my $val=$listArr[$i];
        $iqrSum += $val if(($val >= $iqrLowerBound) and ($val <= $iqrUpperBound));
        $iqrCount++ if(($val >= $iqrLowerBound) and ($val <= $iqrUpperBound));
    }
    $iqrMean=($iqrSum/$iqrCount) if($iqrCount != 0);
    
    my $iqrStdev="NA";
    my $iqrVariance="NA";
    
    if($iqrCount > 1) {
        my $total_iqr_deviation=0;
        for(my $i=0;$i<$listArrSize;$i++) {
            my $val=$listArr[$i];
            
            # only use those values within the new IQR bounds
            next if(($val <= $iqrLowerBound) or ($val >= $iqrUpperBound));
            
            my $deviation=$val-$iqrMean;
            my $sqr_deviation=($deviation**2);
            $total_iqr_deviation += $sqr_deviation;
        }
        $iqrStdev=($total_iqr_deviation/($iqrCount-1));
        $iqrVariance=$iqrStdev;
        $iqrStdev=sqrt($iqrStdev);
    }
        
    $min=$listArr[0]; # get first element (smallest)
    $max=$listArr[-1]; # get last element (largest)
        
    # perform desired trimming
    my $trimmedStart=floor($listArrSize*$trimPC);
    my $trimmedEnd=($listArrSize-$trimmedStart);
    
    #handle small list issue
    $trimmedStart=0 if($trimmedEnd < $trimmedStart);
    $trimmedEnd=($listArrSize) if($trimmedEnd < $trimmedStart);
            
    confess "no elements in the list ($listArrSize) [$trimmedStart - $trimmedEnd]" if(($trimmedEnd-$trimmedStart) < 0);
    
    my @trimmedListArr=();
    for(my $i=$trimmedStart;$i<$trimmedEnd;$i++) {
        my $val=$listArr[$i];
        push(@trimmedListArr,$val);
    }
    
    my $trimmedListArrSize = @trimmedListArr;
    confess "trimmed list is EMPTY!" if($trimmedListArrSize == 0);
    
    my $trimmedMin=0;
    my $trimmedMax=0;
    
    if($trimmedListArrSize > 0) {    
        
        #handle small list issue
        $trimmedStart=0 if($trimmedEnd < $trimmedStart);
        $trimmedEnd=($trimmedListArrSize-1) if($trimmedEnd < $trimmedStart);
                
        confess "no elements in the list ($trimmedListArrSize) [$trimmedStart - $trimmedEnd]!" if(($trimmedEnd-$trimmedStart) < 0);
        
        $trimmedMin=$trimmedListArr[0]; # get first element (smallest)
        $trimmedMax=$trimmedListArr[-1]; # get last element (largest)
        
        $sum=0;
        my $sumLogs=0;
        for(my $i=0;$i<$trimmedListArrSize;$i++) {
            my $val=$trimmedListArr[$i];
            $sum += $val;
                        
            $zeroPC++ if($val == 0);
            
            my $logTmpVal=0;
            $logTmpVal = log($val)/log(2) if($val > 0);
            $skipGeoMean=1 if($val < 0);
            
            $sumLogs += $logTmpVal;
            $count++;
        }
                
        $mean = ($sum / $count);
        $geomean = ($sumLogs / $count);
        $geomean = (2 ** $geomean);
        
        
        my $total_deviation=0;
        for(my $i=0;$i<$trimmedListArrSize;$i++) {
            my $val=$trimmedListArr[$i];
            my $deviation=$val-$mean;
            my $sqr_deviation=$deviation**2;
            $total_deviation=$total_deviation+$sqr_deviation;
        }
        
        if($trimmedListArrSize > 1) {
            $stdev=$total_deviation/($trimmedListArrSize-1);
            $variance=$stdev;
            $stdev=sqrt($stdev);
        }
            
        $binary=1;
    }
    
    #calculate MAD
    my $mad=0;
    
    my @madArr=();
    for(my $i=0;$i<$trimmedListArrSize;$i++) {
        my $val=$trimmedListArr[$i];
        my $absoluteDeviation=abs($val-$median);
        push(@madArr,$absoluteDeviation);        
    }
    @madArr = sort { $a <=> $b } @madArr;
    my $madArrSize=@madArr;
    if(($madArrSize % 2)) {
        my $madIndex = floor($madArrSize * 0.5);
        $mad = $madArr[$madIndex];
    } else {
        my $madIndexLeft = floor($madArrSize * 0.5) - 1;
        my $madIndexRight = $madIndexLeft + 1;        
        my $madLeft = $madArr[$madIndexLeft];
        my $madRight = $madArr[$madIndexRight];
        $mad = (($madLeft + $madRight) / 2);        
    }
    
    $geomean="NA" if($skipGeoMean);
    $zeroPC = round((($zeroPC / $count) * 100),3) if($count != 0);
    
    my %dataHash=();
    $dataHash{ sum }=$sum;
    $dataHash{ count }=$count;
    $dataHash{ zeroPC }=$zeroPC;
    $dataHash{ mean }=$mean;
    $dataHash{ geomean }=$geomean;
    $dataHash{ stdev }=$stdev;
    $dataHash{ variance }=$variance;
    $dataHash{ median }=$median;
    $dataHash{ q1 }=$q1;
    $dataHash{ q3 }=$q3;
    $dataHash{ iqr }=$iqr;
    $dataHash{ iqrMean }=$iqrMean;
    $dataHash{ iqrStdev }=$iqrStdev;
    $dataHash{ iqrVariance }=$iqrVariance;
    $dataHash{ mad }=$mad;
    $dataHash{ min }=$min;
    $dataHash{ max }=$max;
    $dataHash{ trimmedMin }=$trimmedMin;
    $dataHash{ trimmedMax }=$trimmedMax;
    $dataHash{ binary }=$binary;
    
    return(\%dataHash);
    
}

=head2 roundNearest

 Title     : roundNearest
 Usage     : $rounded_num=roundNearest(...)
 Function  : return rounded num
 Returns   : num
 Argument  : num, nearest digit
 Example   : roundNearest(1501,25)=1500
             roundNearest(1524,25)=1520

=cut

sub roundNearest($$) {
    my $num=shift;
    my $nearest=shift;
    
    return(round($num)) if($nearest == 0);
    
    $num=($num/$nearest);
    $num=round($num);
    $num=($num*$nearest);
    
    return($num);
}

=head2 round

 Title     : round
 Usage     : $rounded_num=round(...)
 Function  : return tradionally rounded num
 Returns   : num
 Argument  : num, nearest digit
 Example   : roundNearest(1501.29525,2)=1501.30
             roundNearest(1501.29525,4)=1501.2953

=cut
    
sub round($;$) {
    # required
    my $num=shift;
    # optional
    my $digs_to_cut=0;
    $digs_to_cut = shift if @_;
    
    return($num) if($num eq "NA");
    
    my $roundedNum=$num;
    
    if(($num != 0) and ($digs_to_cut == 0)) {
        $roundedNum = int($num + $num/abs($num*2));
    } else {
        $roundedNum = sprintf("%.".($digs_to_cut)."f", $num) if($num =~ /\d+\.(\d){$digs_to_cut,}/);
    }
    
    return($roundedNum);
}

=head2 getHeaderObject

 Title     : getHeaderObject
 Usage     : $headerObject=getHeaderObject(...)
 Function  : header string to header object
 Returns   : hash
 Argument  : string, flag
 
=cut

sub getHeaderObject($;$) {
    # required
    my $header=shift;
    # optional 
    my $enforceValidHeaders=0;
    $enforceValidHeaders=shift if @_;
    
    my @tmp=();
    my $tmpSize=0;
    
    my ($subName,$assembly,$coords);
    $subName=$assembly=$coords="NA";    
    @tmp=split(/\|/,$header);
    $tmpSize=scalar @tmp;
    ($subName,$assembly,$coords)=split(/\|/,$header) if($tmpSize == 3);    
    confess "header [$header] is not in proper format...!" if(($enforceValidHeaders) and ($tmpSize != 3));
    
    my ($chromosome,$pos);
    $chromosome=$pos="NA";
    @tmp=split(/:/,$coords);
    $tmpSize=scalar @tmp;
    $pos=$tmp[-1];
    $chromosome = $coords;
    $chromosome =~ s/:$pos//;
    confess "coordinates [$coords] are not in proper format...!" if(($enforceValidHeaders) and ($tmpSize != 2));
    
    my ($region);
    $region=$chromosome;
    
    my $primerType="NA";
    if($subName =~ /__/) {
        @tmp=split(/__/,$subName);
        $region=$tmp[0];
    } else {
        @tmp=split(/_/,$subName);
        $tmpSize=scalar @tmp;
        $region=$tmp[1]."_".$tmp[2] if(($tmpSize == 5) and ($tmp[0] eq "5C"));
        $primerType=$tmp[3] if(($tmpSize == 5) and ($tmp[0] eq "5C"));
    }
    
    my ($start,$end);
    $start=$end=0;
    @tmp=split(/-/,$pos);
    $tmpSize=scalar @tmp;
    ($start,$end)=split(/-/,$pos) if($tmpSize == 2);
    confess "position [$pos] is not in proper format...!" if(($enforceValidHeaders) and ($tmpSize != 2));
    
    my $size=(($end-$start)+1); # add to for 1-based positioning
    my $midpoint=(($end+$start)/2);
        
    my %headerObject=();
    $headerObject{ subName }=$subName;
    $headerObject{ primerType }=$primerType;
    $headerObject{ assembly }=$assembly;
    $headerObject{ chromosome }=$chromosome;
    $headerObject{ coords }=$coords;
    $headerObject{ region }=$region;
    $headerObject{ start }=$start;
    $headerObject{ end }=$end;
    $headerObject{ midpoint }=$midpoint;
    $headerObject{ size }=$size;
    
    return(\%headerObject);
    
}

=head2 getInteractionDistance

 Title     : getInteractionDistance
 Usage     : $distance=getInteractionDistance(...)
 Function  : return bp distance between two headers/bins
 Returns   : number
 Argument  : matrixObject hash, headerObject hash 1, headerObject hash 2
 
=cut

sub getInteractionDistance($$$;$$) {
    #required
    my $matrixObject=shift;
    my $headerObject1=shift;
    my $headerObject2=shift;
    #optional
    my $cisApproximateFactor=1;
    $cisApproximateFactor=shift if @_;
    my $logTransform=0;
    $logTransform=shift if @_;

    my $chr_1=$headerObject1->{ chromosome };
    my $region_1=$headerObject1->{ region };
    my $start_1=$headerObject1->{ start };
    my $end_1=$headerObject1->{ end };
    my $midpoint_1=$headerObject1->{ midpoint };
    
    my $chr_2=$headerObject2->{ chromosome };
    my $region_2=$headerObject2->{ region };
    my $start_2=$headerObject2->{ start };
    my $end_2=$headerObject2->{ end };
    my $midpoint_2=$headerObject2->{ midpoint };
    
    return(-1) if(($chr_1 eq "NA") or ($chr_2 eq "NA"));
    return(-1) if(($start_1 >= $end_1) or ($start_2 >= $end_2));
    return(-1) if($region_1 ne $region_2);    
    return(-1) if($chr_1 ne $chr_2);    
    
    my $equalHeaderFlag = $matrixObject->{ equalHeaderFlag };
    my $headerSizing = $matrixObject->{ headerSizing };
    
    # override midpoint, if equalHeaderFlag == 1 (to deal with many many extra distances because of half empty final bin per contig)
    $midpoint_1 = ($start_1+($headerSizing/2)) if($equalHeaderFlag);
    $midpoint_2 = ($start_2+($headerSizing/2)) if($equalHeaderFlag);    
    
    my $dist=-1;
    if($equalHeaderFlag == 0) { # bins do not overlap
        if($midpoint_1 == $midpoint_2) { #self
            $dist = 0;
        } else {
            if($start_1 > $start_2) { 
                $dist = abs($start_1-$end_2);
            } else { 
                $dist = abs($start_2-$end_1);
            }
        }
    } else { # bins do overlap
        $dist = abs($midpoint_1-$midpoint_2);
    }    
    
    #transform dist into approximate dist if necessary
    $dist = round($dist/$cisApproximateFactor) if(($dist != -1) and ($dist != 0)); #do not re-scale if TRANS or SELF
    $dist = log($dist)/log($logTransform) if(($logTransform > 0) and ($dist > 0));
    
    return($dist);
}

=head2 classifyInteraction

 Title     : classifyInteraction
 Usage     : $interactionType=classifyInteraction(...)
 Function  : return type of interaction [USABLE,NULL] 
 Returns   : string
 Argument  : matrixObject hash, includeCis flag, maxDistance, includeTrans flag, headerObject hash 1, headerObject hash 2
 
=cut

sub classifyInteraction($$$$$$$) {
    my $matrixObject=shift;
    my $includeCis=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $includeTrans=shift;
    my $headerObject1=shift;
    my $headerObject2=shift;
    
    # get true interactor distance (non-scaled)
    my $interactionDistance=getInteractionDistance($matrixObject,$headerObject1,$headerObject2);
    
    my $chr_1=$headerObject1->{ chromosome };
    my $region_1=$headerObject1->{ region };
    
    my $chr_2=$headerObject2->{ chromosome };
    my $region_2=$headerObject2->{ region };
    
    return("USABLE") if(($includeTrans) and ($interactionDistance == -1) and ($chr_1 ne $chr_2));
    
    
    if($interactionDistance == -1) { # trans
        return("USABLE") if( $includeTrans );
    } else { #cis
        return("USABLE") if( ( $includeCis ) and (!defined($minDistance)) and (!defined($maxDistance)) ); # handle no distance criteria
        return("NULL") if( ( $includeCis ) and (defined($minDistance)) and ($interactionDistance < $minDistance) ); # handle lower bound
        return("NULL") if( ( $includeCis ) and (defined($maxDistance)) and ($interactionDistance > $maxDistance) ); # handle upper bound
        # remaining = within distance range
        
        return("USABLE");
    }
    
    confess "error classifying interaction [$interactionDistance]";
}

=head2 parseContigs

 Title     : parseContigs
 Usage     : $header2contig,$index2contig,$contig2index,$contigList,$assembly=parseContigs(...)
 Function  : returns set of contig objects from matrix file
 Returns   : hash, hash, hash, hash
 Argument  : input matrix file path, inc2header hash, header2inc hash
 
=cut

sub parseContigs($$$) {
    my $inputMatrix=shift;
    my $inc2header=shift;
    my $header2inc=shift;
    
    my $numYHeaders=keys(%{$header2inc->{ y }});
    my $numXHeaders=keys(%{$header2inc->{ x }});
    my $numTotalHeaders=keys(%{$header2inc->{ xy }});
    
    my $header2contig={};
    my $contig2index={};
    my $index2contig={};
    my $contig2inc={};
    my $header2contiginc={};
    my $contiginc2header={};
    my $contigList={};
    
    my $contigIndex=-1;
    my ($contigInc);
    my $assembly="NA";
    
    my $lastYContig="NA";
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        my $yContig=$yHeaderObject->{ region };
        my $yAssembly=$yHeaderObject->{ assembly };
        $assembly=$yAssembly if($assembly eq "NA");
        
        #confess "non-constant assembly! [$yAssembly:$assembly]" if($yAssembly ne $assembly);
        
        $contigIndex = $contig2index->{$yContig} if(exists($contig2index->{$yContig}));
        $contigIndex = (keys %{$contig2index->{ xy }}) if(!exists($contig2index->{ xy }->{$yContig}));
        
        $contigInc = $header2contiginc->{ xy }->{$yHeader} if(exists($header2contiginc->{ xy }->{$yHeader}));
        $contigInc = $contig2inc->{$yContig}++ if(!exists($header2contiginc->{ xy }->{$yHeader}));
        
        $contig2index->{ y }->{$yContig}=$contigIndex;
        $index2contig->{ y }->{$contigIndex}=$yContig;
        $contig2index->{ xy }->{$yContig}=$contigIndex;
        $index2contig->{ xy }->{$contigIndex}=$yContig;
        
        $header2contig->{ y }->{$yHeader}=$yContig;
        $header2contiginc->{ y }->{$yHeader}=$contigInc;
        $header2contig->{ xy }->{$yHeader}=$yContig;
        $header2contiginc->{ xy }->{$yHeader}=$contigInc;
        
        $contiginc2header->{$yContig}->{$contigInc}=$yHeader;
        $lastYContig=$yContig;
    }
    
    my $lastXContig="NA";
    for(my $x=0;$x<$numXHeaders;$x++) {
        my $xHeader=$inc2header->{ x }->{$x};
        my $xHeaderObject=getHeaderObject($xHeader);
        my $xContig=$xHeaderObject->{ region };
        my $xAssembly=$xHeaderObject->{ assembly };
        $assembly=$xAssembly if($assembly eq "NA");

        #confess "non-constant assembly! [$xAssembly:$assembly]" if($xAssembly ne $assembly);
        
        $contigIndex = $contig2index->{$xContig} if(exists($contig2index->{$xContig}));
        $contigIndex = (keys %{$contig2index->{ xy }}) if(!exists($contig2index->{ xy }->{$xContig}));
        
        $contigInc = $header2contiginc->{ xy }->{$xHeader} if(exists($header2contiginc->{ xy }->{$xHeader}));
        $contigInc = $contig2inc->{$xContig}++ if(!exists($header2contiginc->{ xy }->{$xHeader}));
        
        $contig2index->{ x }->{$xContig}=$contigIndex;
        $index2contig->{ x }->{$contigIndex}=$xContig;
        $contig2index->{ xy }->{$xContig}=$contigIndex;
        $index2contig->{ xy }->{$contigIndex}=$xContig;
        
        $header2contig->{ x }->{$xHeader}=$xContig;
        $header2contiginc->{ x }->{$xHeader}=$contigInc;
        $header2contig->{ xy }->{$xHeader}=$xContig;
        $header2contiginc->{ xy }->{$xHeader}=$contigInc;
        
        $contiginc2header->{$xContig}->{$contigInc}=$xHeader;
        $lastXContig=$xContig;
    }
    
    my $numContigs=keys(%{$contig2index->{ xy }});
    
    for(my $c=0;$c<$numContigs;$c++) {
        my $contig=$index2contig->{ xy }->{$c};
        
        my $nContigHeaders=$contig2inc->{$contig};
        for(my $ci=0;$ci<$nContigHeaders;$ci++) {
            my $contigHeader=$contiginc2header->{$contig}->{$ci};

            my $contigHeaderObject=getHeaderObject($contigHeader);
            my $contigHeaderStart=$contigHeaderObject->{ start };
            my $contigHeaderEnd=$contigHeaderObject->{ end };
            my $contigHeaderChromosome=$contigHeaderObject->{ chromosome };
            my $contigHeaderAssembly=$contigHeaderObject->{ assembly };
            
            $contigList->{$contig}->{ contigStart }=$contigHeaderStart if( (!exists($contigList->{$contig}->{ contigStart })) or ($contigHeaderStart < $contigList->{$contig}->{ contigStart }) );
            $contigList->{$contig}->{ contigEnd }=$contigHeaderEnd if( (!exists($contigList->{$contig}->{ contigEnd })) or ($contigHeaderEnd > $contigList->{$contig}->{ contigEnd }) );
            $contigList->{$contig}->{ contigAssembly }=$contigHeaderAssembly;
            $contigList->{$contig}->{ contigChromosome }=$contigHeaderChromosome;
            $contigList->{$contig}->{ contigLength }=($contigList->{$contig}->{ contigEnd }-$contigList->{$contig}->{ contigStart });
        }
        my $contigStart=$contigList->{$contig}->{ contigStart };
        my $contigEnd=$contigList->{$contig}->{ contigEnd };
        my $contigAssembly=$contigList->{$contig}->{ contigAssembly };
        my $contigChromosome=$contigList->{$contig}->{ contigChromosome };
        my $contigLength=$contigList->{$contig}->{ contigLength };
    }
    
    $assembly=stripAssemblyGroup($assembly);
    
    return($header2contig,$index2contig,$contig2index,$contigList,$assembly);
    
}

=head2 validateMatrixFile

 Title     : validateMatrixFile
 Usage     : validateMatrixFile(...)
 Function  : ensures file path is valid
 Returns   : None
 Argument  : input matrix file path
 
=cut

sub validateMatrixFile($) {
    my $inputMatrix=shift;
        
    return(1) if(($inputMatrix =~ /\.gz$/) and (!(-T($inputMatrix))));
    return(1) if((-T($inputMatrix)) and ($inputMatrix !~ /.png$/) and ($inputMatrix !~ /.gz$/));
    
    return(0);
}

=head2 checkHeaders

 Title     : checkHeaders
 Usage     : $flag=checkHeaders(...)
 Function  : ensure matrix headers are valid
 Returns   : Flag
 Argument  : input matrix file path
 
=cut

sub checkHeaders($) {
    #required
    my $inputMatrix=shift;

    my $headerFlag=1;
    my $lineNum=0;
    my $init=1;
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    
    while(my $line = <IN>) {
        chomp($line);
        next if($line =~ /^#/);
        next if(($line eq "") or ($line =~ m/^#/));
        
        last if($lineNum > 0);
        
        if($lineNum == 0) {
            my @xHeaders=split(/\t/,$line);
            my $xhsize=@xHeaders;
            
            $init=0;
            my %tmpXHeaders=();
            for(my $x=1;$x<$xhsize;$x++) { # skip to left of matrix
                my $xHead=$xHeaders[$x];
                $headerFlag = 0 if($xHead =~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                $headerFlag = 0 if(exists($tmpXHeaders{$xHead})); # if duplicate header, then use row/col # as header
                $tmpXHeaders{$xHead}=1;
            }
            undef %tmpXHeaders;
        }
        $lineNum++;
        #last;
    }
    
    close(IN);
    
    return($headerFlag);
        
}

=head2 parseHeaders

 Title     : parseHeaders
 Usage     : $inc2header,$header2inc=parseHeaders(...)
 Function  : extract headers from input matrix
 Returns   : inc2header hash, header2inc hash
 Argument  : input matrix file path
 
=cut

sub parseHeaders($) {
    #required
    my $inputMatrix=shift;
    
    my $inc2header={};    
    my $header2inc={};
    
    my $noHeaderFlag=0;
    my $headerCornerFlag=0;
    
    my %tmpXHeaders=();
    my %tmpYHeaders=();
    
    my ($lineNum,$numXHeaders,$numYHeaders);
    $lineNum=$numXHeaders=$numYHeaders=0;
    
    # subtract 1 to get rid of header line
    my $numLines=getNumberOfLines($inputMatrix)-1;
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            my @xHeaders=split(/\t/,$line);
            my $xhsize=@xHeaders;
            
            $headerCornerFlag = 1 if(($xHeaders[0] eq "") or ($xHeaders[0] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/)));
            
            for(my $x=0;$x<$xhsize;$x++) {
                my $xHead=$xHeaders[$x];
                $noHeaderFlag = 1 if($xHead =~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                confess "non-unique x headers!" if(exists($tmpYHeaders{$xHead}));
                $tmpXHeaders{$xHead}=1;
            }
            undef %tmpXHeaders;
            
            for(my $x=1;$x<$xhsize;$x++) {
                my $xHead=$xHeaders[$x];
                $xHead="x".$numXHeaders if($noHeaderFlag);
                $header2inc->{ x }->{$xHead}=$numXHeaders;
                $inc2header->{ x }->{$numXHeaders}=$xHead;
                $numXHeaders++;
            }
            
        } else {
            
            my @data=split(/\t/,$line);
            my $dsize=@data;
            my $yHead=$data[0];
            
            # if X is > 10000, assume symmetrical
            if(($numXHeaders > 10000) and ($numLines - ($numXHeaders-1) >= 0) and ($yHead eq $inc2header->{ x }->{0})) {
                close(IN);
                $header2inc->{ y }=$header2inc->{ x };
                $inc2header->{ y }=$inc2header->{ x };
                $header2inc->{ xy }=$header2inc->{ x };
                $inc2header->{ xy }=$inc2header->{ x };
                return($inc2header,$header2inc);
            }
            
            $noHeaderFlag = 1 if($yHead =~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
            confess "non-unique y headers!" if(exists($tmpYHeaders{$yHead}));
            
            $yHead="y".$numYHeaders if($noHeaderFlag);
            $header2inc->{ y }->{$yHead}=$numYHeaders;
            $inc2header->{ y }->{$numYHeaders}=$yHead;
            
            $numYHeaders++;
            
            confess "matrix too large! > 20000 row/cols" if($numYHeaders > 20000);
        }
        $lineNum++;
    }
    
    undef %tmpYHeaders;
    
    my $symmetrical=isSymmetrical($inc2header);
    if($symmetrical) {
        $header2inc->{ xy }=$header2inc->{ y };
        $inc2header->{ xy }=$inc2header->{ y };
        return($inc2header,$header2inc);
    }
    
    # combine and sort the headers
    my %uniqueHeaders=();
    my @headers=();
    my $h=0;
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        my $yHeaderStart=$yHeaderObject->{ start };
        
        my $yHeaderChromosome=header2subMatrix($yHeader,'liteChr');
        my $yHeaderGroup=header2subMatrix($yHeader,'group');
        
        next if(exists($uniqueHeaders{$yHeader}));
        
        $headers[$h]{ group } = $yHeaderGroup;
        $headers[$h]{ chr } = $yHeaderChromosome;
        $headers[$h]{ start } = $yHeaderStart;
        $headers[$h]{ header } = $yHeader;
        $uniqueHeaders{$yHeader}=1;
        $h++;
    }
    for(my $x=0;$x<$numXHeaders;$x++) {
        my $xHeader=$inc2header->{ x }->{$x};
        my $xHeaderObject=getHeaderObject($xHeader);
        my $xHeaderStart=$xHeaderObject->{ start };
        
        my $xHeaderChromosome=header2subMatrix($xHeader,'liteChr');
        my $xHeaderGroup=header2subMatrix($xHeader,'group');
        
        next if(exists($uniqueHeaders{$xHeader}));
        
        $headers[$h]{ group } = $xHeaderGroup;
        $headers[$h]{ chr } = $xHeaderChromosome;
        $headers[$h]{ start } = $xHeaderStart;
        $headers[$h]{ header } = $xHeader;
        $uniqueHeaders{$xHeader}=1;
        $h++;
    }    
    
    @headers = sort { $a->{ header } cmp $b->{ header } } @headers;
    @headers = sort { $a->{ start } <=> $b->{ start } } @headers;
    @headers = sort { $a->{ chr } cmp $b->{ chr } } @headers;
    @headers = sort { $a->{ group } cmp $b->{ group } } @headers;
    
    for(my $i=0;$i<@headers;$i++) {
        my $header=$headers[$i]{ header };
        if( (!exists($header2inc->{ xy }->{$header})) and (!exists($inc2header->{ xy }->{$i})) ) {
            $header2inc->{ xy }->{$header}=$i;
            $inc2header->{ xy }->{$i}=$header;
        }
    }
        
    close(IN);
    
    return($inc2header,$header2inc);
        
}

=head2 updateMatrixObject

 Title     : updateMatrixObject
 Usage     : $matrixObject=updateMatrixObject(...)
 Function  : update matrix object
 Returns   : matrixObject
 Argument  : matrixObject hash
 
=cut

sub updateMatrixObject($) {
    my $matrixObject=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=keys(%{$header2inc->{ y }});
    my $numXHeaders=keys(%{$header2inc->{ x }});
    
    undef $header2inc->{ xy };
    undef $inc2header->{ xy };
    
    $matrixObject->{ numYHeaders }=$numYHeaders;
    $matrixObject->{ numXHeaders }=$numXHeaders;
    $matrixObject->{ numTotalHeaders }="NA";
    
    # combine and sort the headers
    my %uniqueHeaders=();
    my @headers=();
    my $h=0;
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        my $yHeaderChromosome=$yHeaderObject->{ chromosome };
        my $yHeaderStart=$yHeaderObject->{ start };
        
        next if(exists($uniqueHeaders{$yHeader}));
        
        $headers[$h]{ chr } = $yHeaderChromosome;
        $headers[$h]{ start } = $yHeaderStart;
        $headers[$h]{ header } = $yHeader;
        $uniqueHeaders{$yHeader}=1;
        $h++;
    }
    for(my $x=0;$x<$numXHeaders;$x++) {
        my $xHeader=$inc2header->{ x }->{$x};
        my $xHeaderObject=getHeaderObject($xHeader);
        my $xHeaderChromosome=$xHeaderObject->{ chromosome };
        my $xHeaderStart=$xHeaderObject->{ start };
        
        next if(exists($uniqueHeaders{$xHeader}));
        
        $headers[$h]{ chr } = $xHeaderChromosome;
        $headers[$h]{ start } = $xHeaderStart;
        $headers[$h]{ header } = $xHeader;
        $uniqueHeaders{$xHeader}=1;
        $h++;
    }    
    
    @headers = sort { $a->{ header } cmp $b->{ header } } @headers;
    @headers = sort { $a->{ start } <=> $b->{ start } } @headers;
    @headers = sort { $a->{ chr } cmp $b->{ chr } } @headers;
    
    for(my $i=0;$i<@headers;$i++) {
        my $header=$headers[$i]{ header };
        if( (!exists($header2inc->{ xy }->{$header})) and (!exists($inc2header->{ xy }->{$i})) ) {
            $header2inc->{ xy }->{$header}=$i;
            $inc2header->{ xy }->{$i}=$header;
        }
    }
    
    my $numTotalHeaders=keys(%{$header2inc->{ xy }});
    
    $matrixObject->{ inc2header }=$inc2header;
    $matrixObject->{ header2inc }=$header2inc;
    
    my $yMaxHeaderLength=getMaxHeaderLength($inc2header->{ y });
    my $xMaxHeaderLength=getMaxHeaderLength($inc2header->{ x });
    
    my $symmetrical=isSymmetrical($inc2header);
    
    # calculate number of interactions
    my $numInteractions=($numYHeaders*$numXHeaders);
    $numInteractions=((($numTotalHeaders*$numTotalHeaders)-$numTotalHeaders)/2) if($symmetrical);
    
    $matrixObject->{ numInteractions }=$numInteractions;
    $matrixObject->{ xHeaderLength }=$xMaxHeaderLength;
    $matrixObject->{ yHeaderLength }=$yMaxHeaderLength;
    $matrixObject->{ symmetrical }=$symmetrical;
    
    return($matrixObject);
}
    
=head2 getNArows

 Title     : getNArows
 Usage     : $NA_rowcols=getNArows(...)
 Function  : find NA row / cols
 Returns   : NA_rowcols hash
 Argument  : input matrix file path
 
=cut

sub getNARows($) {
    #required
    my $inputMatrix=shift;
    
    my $lineNum=0;
    my %NA_headers=();
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum > 0) { # skip x headers
            my @data=split(/\t/,$line);
            my $dsize=@data;
            my $yHead=$data[0];
            
            my $naCount=0;
            for(my $d=1;$d<$dsize;$d++) {
                my $score=$data[$d];
                last if($score ne "NA");
                $naCount++;
            }
            
            $NA_headers{$yHead}=1 if($naCount == ($dsize-1));
        }
        $lineNum++;
    }
    close(IN);
    
    return(\%NA_headers);
}

=head2 logTransformMatrix

 Title     : logTransformMatrix
 Usage     : $matrix=logTransformMatrix(...)
 Function  : log transform interaction matrix
 Returns   : 2D hash
 Argument  : matrix 2D hash, matrixObject hash, logTransform number
 
=cut

sub logTransformMatrix($$$) {
    my $matrix=shift;
    my $matrixObject=shift;
    my $logTransform=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $cScore=$missingValue;
            $cScore = $matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            if($logTransform > 0) {
                my $tmp_cScore = "NA";
                $tmp_cScore = (log($cScore)/log($logTransform)) if($cScore > 0);
                $cScore=$tmp_cScore;
            }
            
            $matrix->{$y}->{$x}=$cScore;
        }
    }
    
    return($matrix);
}

=head2 getData

 Title     : getData
 Usage     : $matrix,$matrixObject=getData(...)
 Function  : matrix file to 2D hash, matrixObject
 Returns   : matrix 2D hash, matrixObject hahs
 Argument  : input matrix file path, matrixObject hash
 
=cut

sub getData($$;$$$$$$) {
    # required
    my $inputMatrix=shift;
    my $matrixObject=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    my $minDistance = undef;
    $minDistance=shift if @_;
    my $maxDistance = undef;
    $maxDistance=shift if @_;
    my $excludeCis=0;
    $excludeCis=shift if @_;
    my $excludeTrans=0;
    $excludeTrans=shift if @_;
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    print STDERR "loading matrix data ...\n" if($verbose);
    
    my $header2inc=$matrixObject->{ header2inc };
    
    my $subsetMode=0;
    $subsetMode = 1 if((defined($maxDistance)) or ($excludeTrans) or ($excludeCis));
    $matrixObject->{ missingValue }="NA" if($subsetMode);
    
    my $symmetrical=$matrixObject->{ symmetrical };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $missingValue=$matrixObject->{ missingValue };
    my $headerFlag=$matrixObject->{ headerFlag };
    
    my (%matrix);
    
    my $lineNum=0;
    my @xHeaders=();
    
    print STDERR "\tgetData\n" if($verbose);
    
    my $nLines = getNumberOfLines($inputMatrix)-1;
    my $progressBucketSize=ceil($nLines / 1000);
    my $pcComplete=0;
    
    my $maxNonZeros = (15000*15000);
    my $nDataPoints=0;
    
    my %headerObjects=();
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            @xHeaders=split(/\t/,$line);
        } else {
            my @data=split(/\t/,$line);
            my $dsize=@data;
            
            my $yHeader=$data[0];
            my $yHeaderObject={};
            $yHeaderObject=getHeaderObject($yHeader) if(($subsetMode) and (!exists($headerObjects{$yHeader})));
            $yHeaderObject=$headerObjects{$yHeader} if(($subsetMode) and (exists($headerObjects{$yHeader})));
            $headerObjects{$yHeader}=$yHeaderObject if(($subsetMode) and (!exists($headerObjects{$yHeader})));
            
            my $yIndex=-1;
            $yIndex = $header2inc->{ y }->{$yHeader} if(defined($header2inc->{ y }->{$yHeader}));
            
            next if($yIndex == -1);
            
            my $indexStart=1;
            my $indexEnd=$dsize;
            
            for(my $i=$indexStart;$i<$indexEnd;$i++) {
                my $cScore=$data[$i];
                
                # skip if cScore is not a valid number
                $cScore = "NA" if($cScore !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                
                # sparse matrix logic - do not store 0/nan (depending on quantity)
                next if( ($cScore eq $missingValue) or (($cScore ne "NA") and ($missingValue ne "NA") and ($cScore == $missingValue)) );
                
                # truncate numbers to minimal digits
                $cScore = sprintf "%.".$sigDigits."f", $cScore if(($cScore ne "NA") and ($cScore !~ /^[+-]?\d+$/));
                
                my $xHeader=$xHeaders[$i];
                my $xHeaderObject={};
                $xHeaderObject=getHeaderObject($xHeader) if(($subsetMode) and (!exists($headerObjects{$xHeader})));
                $xHeaderObject=$headerObjects{$xHeader} if(($subsetMode) and (exists($headerObjects{$xHeader})));
                $headerObjects{$xHeader}=$xHeaderObject if(($subsetMode) and (!exists($headerObjects{$xHeader})));
                
                my $xIndex=-1;
                $xIndex = $header2inc->{ x }->{$xHeader} if(defined($header2inc->{ x }->{$xHeader}));
                next if($xIndex == -1);
                                
                if($subsetMode) {
                    my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
                    next if(($interactionDistance == -1) and ($excludeTrans));
                    next if(($interactionDistance != -1) and ($excludeCis));
                    next if((defined($maxDistance)) and ($interactionDistance != -1) and ($interactionDistance >= $maxDistance));
                }
                
                # ensure symmetrical data
                if(($symmetrical) and (exists($matrix{$xIndex}{$yIndex}))) {
                    confess "data is not symmetrical (x,y) ($xIndex,$yIndex) [$cScore != $matrix{$xIndex}{$yIndex}], yet headers are" if(($matrix{$xIndex}{$yIndex} ne $cScore) and ($subsetMode == 0));
                }
                
                confess "error - matrix is too large ($nDataPoints > $maxNonZeros)!" if($nDataPoints > $maxNonZeros);
                
                $matrix{$yIndex}{$xIndex}=$cScore;
                $nDataPoints++
                
            }
        }
                
        $pcComplete = 100 if($lineNum == ($nLines));
        print STDERR "\e[A" if(($verbose) and ($lineNum != 0));
        printf STDERR "\t%.2f%% complete ($lineNum/".($nLines).")...\n", $pcComplete if($verbose);
        $pcComplete = round((($lineNum/$nLines)*100),2) if($nLines != 0);
        $lineNum++;
    }
    close(IN);
    
    print STDERR "\tloaded ".commify($nDataPoints)." datapoints\n" if($verbose);
    
    confess "no data available!" if($nDataPoints == 0);
    
    return(\%matrix,$matrixObject);
}

=head2 readLoessFile

 Title     : readLoessFile
 Usage     : $loess=readLoessFile(...)
 Function  : read loess file into loess hash
 Returns   : loess hash
 Argument  : loes file path
 
=cut

sub readLoessFile($;$) {
    # required
    my $loessFile=shift;
    # optional
    my $verbose=0;
    $verbose = shift if @_;
    
    confess "loess file does not exist [$loessFile]" if(!(-e $loessFile));
        
    my %loess=();
    
    my $lineNum=0;
    my %header2index=();
    
    print STDERR "\trecovering loessFile\n" if($verbose);
    
    open(IN,inputWrapper($loessFile)) or confess "Could not open file [$loessFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            my @headers=split(/\t/,$line);
            for(my $i=0;$i<@headers;$i++) {
                my $header=$headers[$i];
                $header2index{$header}=$i;
            }
            $lineNum++;
            next;
        }
        
        my @tmp=split(/\t/,$line);
        
        my $interactionDistance = $tmp[ $header2index{ interactionDistance } ];
        my $loessValue = $tmp[ $header2index{ loessExpectedValue } ];
        my $loessStdev = $tmp[ $header2index{ loessExpectedStdev } ];
        
        $loess{$interactionDistance}{ loess }=$loessValue;
        $loess{$interactionDistance}{ stdev }=$loessStdev;
        
        $lineNum++;
        
    }
    
    close(IN);
    
    print STDERR "\n" if($verbose);
    
    return(\%loess);
}    

=head2 validateLoessObject

 Title     : validateLoessObject
 Usage     : validateLoessObject(...)
 Function  : validate loess object file
 Returns   : bool
 Argument  : loess object file 
 
=cut

sub validateLoessObject($) {
    # required
    my $loessObjectFile=shift;
    
    return(0) if(!(-e $loessObjectFile));
    
    my $last="";
    open(IN,inputWrapper($loessObjectFile));
    while (<IN>) { $last = $_ }
    chomp($last);
   
    return(1) if($last eq "## done");
    
    return(0);
}

=head2 calculateLoess

 Title     : calculateLoess
 Usage     : $loess=calculateLoess(...)
 Function  : calculate expected [loess] on matrix
 Returns   : loess hash
 Argument  : matrixObject hash, inputData arr, loessFile, loess alpha value [0-1], disableLoessIQRFilter flag, skipZero flag
 
=cut

sub calculateLoess($$$$$$$;$) {
    # required
    my $matrixObject=shift;
    my $inputDataRef=shift;
    my @inputData = @$inputDataRef;
    my $totalDataSize=@inputData;
    my $loessFile=shift;
    my $loessObjectFile=shift;
    my $loessAlpha=shift;
    my $disableLoessIQRFilter=shift;
    my $skipZero=shift;
    # optional 
    my $debug=0;
    $debug=shift if @_;

    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numInteractions=$matrixObject->{ numInteractions };    
    my $verbose=$matrixObject->{ verbose };
    
    
    print STDERR "calculating expected ...\n" if(($verbose) and (validateLoessObject($loessObjectFile)));
    return(readLoessFile($loessObjectFile,$verbose)) if(validateLoessObject($loessObjectFile));
    
    # ensure available cis data
    return({}) if($totalDataSize == 0);
    
    print STDERR "calculating cis-expected ...\n" if($verbose);
        
    my @tmp=split(/\//,$loessFile);
    my $loessFileName=$tmp[-1];
    my $writeLoessFileFlag=1;
    $writeLoessFileFlag=0 if(($loessFile eq "NA") or ($loessFile eq ""));
    
    my $loessSearchSpace=floor($totalDataSize*$loessAlpha);
    my $bucket=1;
    $bucket = floor($totalDataSize/1000) if($totalDataSize > 1000);
    
    my $pcComplete=0;
    my $loessAlphaChangeAttempts=0;
    
    while(($loessSearchSpace <= 3) and ($loessAlphaChangeAttempts < 100)) {
        my $newLoessAlpha = ($loessAlpha * 2);
        print STDERR "\tERROR-9 : loessAlpha is too small ($loessSearchSpace / $totalDataSize)\n\tIncreasing loessAlpha to $loessAlpha -> $newLoessAlpha and retrying...\n";
        $loessSearchSpace=ceil($totalDataSize*$newLoessAlpha);
        $loessAlpha=$newLoessAlpha;
        $loessAlphaChangeAttempts++;
    }
    while(($loessSearchSpace > 100000) and ($loessAlphaChangeAttempts < 100)) {
        my $newLoessAlpha = ($loessAlpha * 0.75);
        $loessSearchSpace=ceil($totalDataSize*$newLoessAlpha);
        print STDERR "\tERROR-9 : loessAlpha is too large ($loessSearchSpace / $totalDataSize)\n\tcannot touch more than 10000 data points\n\tShrinking loessAlpha to $loessAlpha -> $newLoessAlpha and retrying...\n";
        $loessAlpha=$newLoessAlpha;
        $loessAlphaChangeAttempts++;
    }    
    
    print STDERR "\tcalculateLoess\n" if($verbose);
    print STDERR "\tperforming loess over ($loessSearchSpace/$totalDataSize) datapoints...\n" if($verbose);
    
    open(LOESS,outputWrapper($loessFile)) or confess "Could not open file [$loessFile] - $!" if($writeLoessFileFlag);
    print LOESS "yHeaderName\txHeaderName\tinteractionDistance\trealInteractionDistance\tobservedSignal\tloessExpectedValue\tloessExpectedStdev\tzScore\n" if($writeLoessFileFlag);
    
    open(COLLAPSED,outputWrapper($loessObjectFile)) or confess "Could not open file [$loessObjectFile] - $!";
    print COLLAPSED "interactionDistance\trealInteractionDistance\tloessExpectedValue\tloessExpectedStdev\n";
    
    my $debugFile=$loessFile.".debug";
    open(DEBUG,outputWrapper($debugFile)) or confess "Could not open file [$debugFile] - $!" if($debug);
    print DEBUG "anchorApproximateX\tanchorX\t[loessLowerBound - loessUpperBound]\tnDataPoints_1\tzeroPC_1\t[loessLowerBound - loessUpperBound]\tnDataPoints_2\tzeroPC_2\tloessValue\tloessStdev\n" if($debug);
    
    my %loess=();
    
    for(my $i=0;$i<$totalDataSize;$i++) {
        
        if(($i % $bucket) == 0) { 
            if($pcComplete <= 100) {
                print STDERR "\e[A"  if(($verbose) and ($i != 0));
                printf STDERR "\t%.2f%% complete ($i/$totalDataSize)...\n",$pcComplete if($verbose);
                $pcComplete = round((($i/$totalDataSize)*100),2);
            }
        }
        
        my $anchorApproximateX=$inputData[$i][0];
        my $anchorY=$inputData[$i][1];
        my $anchorX=$inputData[$i][3];
        
        my $anchorInteraction=$inputData[$i][2];
        my ($yHeader,$xHeader)=split(/___/,$anchorInteraction);
        
        # only calculate for NEW distances (x valuess)
        if(!exists($loess{$anchorApproximateX})) {
                
            # calculate N% closest data points
            my $leftWall=($i-$loessSearchSpace);
            my $rightWall=($i+$loessSearchSpace);
                        
            my $leftOffset=0;
            $leftOffset = abs($leftWall) if($leftWall < 0);
            my $rightOffset=0;
            $rightOffset = ($rightWall-$totalDataSize+1) if($rightWall >= $totalDataSize);
                        
            $leftWall -= $rightOffset;
            $rightWall += $leftOffset;
                        
            $leftWall = 0 if($leftWall < 0);
            $rightWall=($totalDataSize-1) if($rightWall >= $totalDataSize);
            
            my $wallSpan=($rightWall-$leftWall)+1;
            my $expectedWallSpan=($loessSearchSpace+$loessSearchSpace+1);
                        
            confess "incorrect loess search space ($wallSpan vs $expectedWallSpan)" if($wallSpan != $expectedWallSpan);
            
            # if left wall does not end on a distance break, increase left wall to include all itx of same dist.
            while((($leftWall-1) >= 0) and ($inputData[$leftWall][0] == $inputData[$leftWall-1][0])) {
                $leftWall--;
            }
            
            # if right wall does not end on a distance break, increase right wall to include all itx of same dist.
            while((($rightWall+1) < $totalDataSize) and ($inputData[$rightWall][0] == $inputData[$rightWall+1][0])) {
                $rightWall++;
            }
            
            # fill the tmp array with N% closest data points
            my @scores=();
            my $scoreSize=0;
            for(my $inputArrIndex=$leftWall;$inputArrIndex<=$rightWall;$inputArrIndex++) {
                my $x=$inputData[$inputArrIndex][3];
                my $approximateX=$inputData[$inputArrIndex][0];
                my $y=$inputData[$inputArrIndex][1];
                my $tmpKey=$inputData[$inputArrIndex][2];
                
                my $dval=dist($x,$anchorX);
                
                $scores[$scoreSize]{ dist2ref }=$dval;
                $scores[$scoreSize]{ inputArrIndex }=$inputArrIndex;
                $scores[$scoreSize]{ x }=$x;
                $scores[$scoreSize]{ y }=$y;
                $scores[$scoreSize]{ key }=$tmpKey;
                $scoreSize++;
            }

            # sort by inputArrIndex then dist2ref
            @scores = sort { $a->{ inputArrIndex } <=> $b->{ inputArrIndex } } @scores;
            @scores = sort { $a->{ dist2ref } <=> $b->{ dist2ref } } @scores;
            
            my $tmpLoessSearchSpace=$loessSearchSpace;
                        
            # if search window ends on a non distance break - increase search window
            while(($tmpLoessSearchSpace < ($scoreSize-1)) and ($scores[$tmpLoessSearchSpace]{ dist2ref } == $scores[$tmpLoessSearchSpace+1]{ dist2ref })) {
                $tmpLoessSearchSpace++;
            }
            
            my ($loessValue,$loessStdev,$loessLowerBound,$loessUpperBound,$R2,$nDataPoints_1,$zeroPC_1,$nDataPoints_2,$zeroPC_2);
            $loessValue=$loessStdev=$loessLowerBound=$loessUpperBound=$R2=$nDataPoints_1=$zeroPC_1=$nDataPoints_2=$zeroPC_2="NA";
                        
            # calculate initial loess
            print DEBUG "$anchorApproximateX\t$anchorX" if($debug);
            ($loessValue,$loessStdev,$loessLowerBound,$loessUpperBound,$R2,$nDataPoints_1,$zeroPC_1)=doLoessIteration(\@scores,$tmpLoessSearchSpace,$skipZero,$loessLowerBound,$loessUpperBound,1,$anchorX,0,$verbose);
            # $perform second pass loess (IQR ROBUST)
            print DEBUG "\t[$loessLowerBound - $loessUpperBound]\t$nDataPoints_1\t$zeroPC_1" if($debug);
            ($loessValue,$loessStdev,$loessLowerBound,$loessUpperBound,$R2,$nDataPoints_2,$zeroPC_2)=doLoessIteration(\@scores,$tmpLoessSearchSpace,$skipZero,$loessLowerBound,$loessUpperBound,$disableLoessIQRFilter,$anchorX,1,$verbose) if($disableLoessIQRFilter == 0);
            print DEBUG "\t[$loessLowerBound - $loessUpperBound]\t$nDataPoints_2\t$zeroPC_2" if($debug);
            print DEBUG "\t$loessValue\t$loessStdev\n" if($debug);
            
            $loess{$anchorApproximateX}{ loess }=$loessValue;
            $loess{$anchorApproximateX}{ stdev }=$loessStdev;
            
            next if($writeLoessFileFlag == 0);
            
            my $zScore = "NA";
            $zScore = (($anchorY-$loessValue)/$loessStdev) if(($loessValue ne "NA") and ($loessStdev != 0));
            
            print LOESS "$yHeader\t$xHeader\t$anchorApproximateX\t$anchorX\t$anchorY\t$loessValue\t$loessStdev\t$zScore\n" if($writeLoessFileFlag);
            
            print COLLAPSED "$anchorApproximateX\t$anchorX\t$loessValue\t$loessStdev\n";
            
                 
        } else {
            
            next if($writeLoessFileFlag == 0);
            
            my $loessValue=$loess{$anchorApproximateX}{ loess };
            my $loessStdev=$loess{$anchorApproximateX}{ stdev };
                        
            my $zScore = "NA";
            $zScore = (($anchorY-$loessValue)/$loessStdev) if(($loessValue ne "NA") and ($loessStdev != 0));
            
            print LOESS "$yHeader\t$xHeader\t$anchorApproximateX\t$anchorX\t$anchorY\t$loessValue\t$loessStdev\t$zScore\n" if($writeLoessFileFlag);
        }
    }

    print STDERR "\e[A" if($verbose);
    printf STDERR "\t%.2f%% complete ($totalDataSize/$totalDataSize)...\n",100 if($verbose);
    print STDERR "\n" if($verbose);
    
    print COLLAPSED "## done\n";
    
    close(COLLAPSED);
    close(LOESS) if($writeLoessFileFlag);
    close(DEBUG) if($debug);
    
    return(\%loess);
    
}

=head2 doLoessIteration

 Title     : doLoessIteration
 Usage     : $loessValue,$loessStdev,$loessLowerBound,$loessUpperBound,$R2,$nDataPoints,$tmpScoresZeroPC=doLoessIteration(...)
 Function  : perform a single loess iteration
 Returns   : $loessValue,$loessStdev,$loessLowerBound,$loessUpperBound,$R2,$nDataPoints,$tmpScoresZeroPC
 Argument  : 
 
=cut

sub doLoessIteration($$$$$$$$;$) {
    my $scores=shift;
    my $tmpLoessSearchSpace=shift;
    my $skipZero=shift;
    my $loessLowerBound=shift || "NA";
    my $loessUpperBound=shift || "NA";
    my $disableLoessIQRFilter=shift;
    my $anchorX=shift;
    my $iterationNumber=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    
    my ($sumWeights,$sumWeightedX,$sumWeightedX2,$sumWeightedY,$sumWeightedXY,$sumSquares,$tmpLoessWindowSize);
    $sumWeights=$sumWeightedX=$sumWeightedX2=$sumWeightedY=$sumWeightedXY=$sumSquares=$tmpLoessWindowSize=0;
    
    my $positiveFlag = 1;
    
    my $maxDist=$scores->[$tmpLoessSearchSpace]->{ dist2ref };
    $maxDist = 1 if($maxDist==0);
        
    my @tmpScores=();
    my %tmpScoresHash=();
    
    if(($disableLoessIQRFilter == 0) and (($loessLowerBound eq "NA") and ($loessUpperBound eq "NA"))) {
        print STDERR "\e[A" if($verbose);
        print STDERR "\tWARNING x=$anchorX (missing IQR bounds [$loessLowerBound - $loessUpperBound])\n\n" if($verbose);
    }
    
    for(my $s=0;$s<=$tmpLoessSearchSpace;$s++) {
        my $dist2ref=$scores->[$s]->{ dist2ref };
        my $inputArrIndex=$scores->[$s]->{ index };
        my $x=$scores->[$s]->{ x };
        my $y=$scores->[$s]->{ y };
                                
        # skip if score is 0 and ignoreZero is selected
        next if(($skipZero) and ($y == 0));
        
        # skip if score is beyond robust limits (outlier!)
        next if(($disableLoessIQRFilter == 0) and (($loessLowerBound ne "NA") and ($loessUpperBound ne "NA")) and (($y < $loessLowerBound) or ($y > $loessUpperBound)));
            
        my $scaledDist=($dist2ref/$maxDist);        
        my $weight=tri_cube($scaledDist);
        my $adjustedScore = ($weight*$y);

        next if($weight == 0);

        push(@tmpScores,$y);
        $tmpScoresHash{$tmpLoessWindowSize}{ y } = $y;
        $tmpScoresHash{$tmpLoessWindowSize}{ x } = $x;
        $tmpScoresHash{$tmpLoessWindowSize}{ w } = $weight;
        $tmpScoresHash{$tmpLoessWindowSize}{ wy } = $adjustedScore;
        
        $positiveFlag=0 if($y < 0);
        
        $sumWeights += $weight;
        $sumWeightedX += ($x*$weight);
        $sumWeightedX2 += (($x**2)*$weight);
        $sumWeightedY += ($y*$weight);
        $sumWeightedXY += (($x*$y)*$weight);
        $tmpLoessWindowSize++;
            
    }
    
    # set loess to NULL if not enough datapoints @ distance
    my $nDataPoints=scalar(@tmpScores);
    return("NA","NA","NA","NA","NA",$nDataPoints,"NA") if($nDataPoints <= 1);
    
    # calculate tmpScore distribution
    my $tmpScoresArrStats=listStats(\@tmpScores);
    my $tmpScoresQ1=$tmpScoresArrStats->{ q1 };
    my $tmpScoresMedian=$tmpScoresArrStats->{ median };
    my $tmpScoresQ3=$tmpScoresArrStats->{ q3 };
    my $tmpScoresIQR=$tmpScoresArrStats->{ iqr };
    my $tmpScoresZeroPC=$tmpScoresArrStats->{ zeroPC };
    # use IQR as robust statistic
    $loessLowerBound = $tmpScoresQ1-(1.5*$tmpScoresIQR);
    $loessUpperBound = $tmpScoresQ3+(1.5*$tmpScoresIQR);
    
    # warn user if percent of un-filtered 0s is > 75%
    if($tmpScoresZeroPC > 75) {
        print STDERR "\e[A" if($verbose);
        print STDERR "\tWARNING x=$anchorX (".$tmpScoresZeroPC."% >75%) - cannot accurately calculate LOWESS expected\n\n" if($verbose);
        return("NA","NA",$loessLowerBound,$loessUpperBound,"NA",$nDataPoints,$tmpScoresZeroPC);
    }
    
    my ($slope,$intercept,$denom);
    $slope=$intercept=$denom=0;
    
    # calculate LOESS denom
    $denom=(($sumWeights*$sumWeightedX2) - ($sumWeightedX**2));
    
    # perform the linear regression
    my $loessValue="NA";
    if($denom == 0) { 
        $loessValue=($sumWeightedY/$tmpLoessWindowSize) if($tmpLoessWindowSize != 0);
        print STDERR "\e[A" if($verbose);
        print STDERR "\tWARNING x=$anchorX (non-distinct X) [".@tmpScores." datapoints] - cannot perform linear regression\n\n" if($verbose);
    } else { 
        $slope=((($sumWeights * $sumWeightedXY) - ($sumWeightedX * $sumWeightedY)) / $denom);
        $intercept=((($sumWeightedX2 * $sumWeightedY) - ($sumWeightedX * $sumWeightedXY)) / $denom);
        $loessValue=(($slope * $anchorX)+$intercept);
    }
        
    # override loess to 0 - if all values were >0 and only the linear fit drops negative
    if(($loessValue < 0) and ($positiveFlag)) {
        print STDERR "\e[A" if($verbose);
        print STDERR "\tWARNING x=$anchorX (linear fit yielded negative value) - overriding LOWESS expected to ($loessValue -> NA)\n\n" if($verbose);
        return("NA","NA",$loessLowerBound,$loessUpperBound,"NA",$nDataPoints,$tmpScoresZeroPC);
    }

    # calculate the LOESS STDEV now
    $sumWeights=$sumWeightedX=$sumWeightedX2=$sumWeightedY=$sumWeightedXY=$sumSquares=$tmpLoessWindowSize=0;
    for(my $s=0;$s<=$tmpLoessSearchSpace;$s++) {
        my $dist2ref=$scores->[$s]->{ dist2ref };
        my $inputArrIndex=$scores->[$s]->{ index };
        my $x=$scores->[$s]->{ x };
        my $y=$scores->[$s]->{ y };
                        
        # skip if score is 0 and ignoreZero is selected
        next if(($skipZero) and ($y == 0));
        
        # skip if score is beyond robust limits (outlier!)
        next if(($disableLoessIQRFilter == 0) and (($y < $loessLowerBound) or ($y > $loessUpperBound)));
        
        my $scaledDist=($dist2ref/$maxDist);
        my $weight=tri_cube($scaledDist);
        my $adjustedScore = ($weight*$y);
        
        next if($weight == 0);
        
        $sumWeights += $weight;
        my $square=$weight*(($y-$loessValue)**2);
        $sumSquares += $square;
        $tmpLoessWindowSize++;
    }
    
    # calculate the loess stdev denom
    my $loessStdevDenom=0;
    $loessStdevDenom=((($tmpLoessWindowSize-1)*($sumWeights))/$tmpLoessWindowSize) if( ($tmpLoessWindowSize != 0) and (((($tmpLoessWindowSize-1)*($sumWeights))/$tmpLoessWindowSize) != 0) );
    my $loessStdev="NA";
    $loessStdev=sqrt($sumSquares/$loessStdevDenom) if($loessStdevDenom != 0);

    my $R2=calculateR2(\%tmpScoresHash,$slope,$intercept,$denom);

    # return loess/bounds
    return($loessValue,$loessStdev,$loessLowerBound,$loessUpperBound,$R2,$nDataPoints,$tmpScoresZeroPC);
}

=head2 calculateR2

 Title     : calculateR2
 Usage     : $R2=calculateR2(...)
 Function  : calculate R^2 for loess scatter subset
 Returns   : R^2 value
 Argument  : 2D hash, slope vlaue, intercept value, denom value
 
=cut

sub calculateR2($$$$) {
    my $tmpScores=shift;
    my $slope=shift;
    my $intercept=shift;
    my $denom=shift;
    
    my $nScores=keys %{$tmpScores};
    return("NA") if($nScores == 0);
    
    my $totY=0;    
    for(my $i_tmp=0;$i_tmp<$nScores;$i_tmp++) {
        my $y=$tmpScores->{$i_tmp}->{ wy };
        $totY += $y;
    }
    
    my $meanY=($totY/$nScores);
    
    my $SStot=0;
    my $SSreg=0;
    my $SSres=0;
    for(my $i_tmp=0;$i_tmp<$nScores;$i_tmp++) {
        my $y=$tmpScores->{$i_tmp}->{ y };
        my $x=$tmpScores->{$i_tmp}->{ x };
        my $w=$tmpScores->{$i_tmp}->{ w };
        
        my $f=$meanY;
        $f=(($slope * $x)+$intercept) if($denom != 0);
        
        $SStot += (($y-$meanY)**2)*$w;
        $SSreg += (($f-$meanY)**2)*$w;
        $SSres += (($y-$f)**2)*$w;    
    }
    
    my $R2 = 0;
    $R2 = (1 - ($SSres/$SStot)) if($SStot != 0);
    
    return($R2);
}

=head2 tri_cube

 Title     : tri_cube
 Usage     : $weight=tri_cube(...)
 Function  : tri cubic function for loess weigh
 Returns   : weight
 Argument  : distance
 
=cut

sub tri_cube($) {
    my $distance=shift;
    $distance=abs($distance);
    
    my $weight=0;
    
    if($distance >= 1) { 
        $weight=0; 
    } else {
        $weight = ((1-($distance**3))**3);
    }
    
    return($weight);
}

=head2 dist

 Title     : dist
 Usage     : $dist=dist(...)
 Function  : calculate abs(dist) between two points
 Returns   : dist
 Argument  : value 1, value 2
 
=cut

sub dist($$) {
    my $p1=shift;
    my $p2=shift;
    my $dist=abs($p2-$p1);
    return $dist;
}

=head2 calculateTransExpected

 Title     : calculateTransExpected
 Usage     : $loess=calculateTransExpected(...)
 Function  : calculate trans expected - embed into loess hash
 Returns   : loess hash
 Argument  : input data hash, skipZero flag, loess hash
 
=cut

sub calculateTransExpected($$$$;$) {
    #required
    my $inputData=shift;
    my $skipZero=shift;
    my $loess=shift;
    my $loessObjectFile=shift;
    #optional
    my $verbose=1;
    $verbose=shift if @_;

    my $inputDataSize=scalar(@{$inputData});
    
    # return if no supplied trans data
    return($loess) if($inputDataSize == 0);

    print STDERR "calculating trans-expected ...\n" if($verbose);
    
    my @tmpList=();
    for(my $i=0;$i<$inputDataSize;$i++) {
        my $value=$inputData->[$i]->[1];
        next if(($value == 0) and ($skipZero));
        push(@tmpList,$value);
    }
    
    # perform percentile trimming
    my $tmpArrStats=listStats(\@tmpList,0.10);
    my $tmpMean=$tmpArrStats->{ mean };
    my $tmpStdev=$tmpArrStats->{ stdev };
    my $tmpQ1=$tmpArrStats->{ q1 };
    my $tmpMedian=$tmpArrStats->{ median };
    my $tmpQ3=$tmpArrStats->{ q3 };
    my $tmpIQR=$tmpArrStats->{ iqr };
    my $tmpMin=$tmpArrStats->{ min };
    my $tmpMax=$tmpArrStats->{ max };
    my $tmpIQRMean=$tmpArrStats->{ iqrMean };
    my $tmpIQRStdev=$tmpArrStats->{ iqrStdev };
    
    if($verbose) {
        print STDERR "\tTRANS\tmin\t$tmpMin\n";
        print STDERR "\tTRANS\tmax\t$tmpMax\n";
        print STDERR "\tTRANS\tmean\t$tmpMean\n";
        print STDERR "\tTRANS\tstdev\t$tmpStdev\n";
        print STDERR "\tTRANS\tQ1\t$tmpQ1\n";
        print STDERR "\tTRANS\tQ2 (median)\t$tmpMedian\n";
        print STDERR "\tTRANS\tQ3median\t$tmpQ3\n";
        print STDERR "\tTRANS\tIQR\t$tmpIQR\n";
        print STDERR "\tTRANS\tIQRMean\t$tmpIQRMean\n";
        print STDERR "\tTRANS\tIQRStdev\t$tmpIQRStdev\n";
        print STDERR "\n" if($verbose);
    }
    
    $loess->{-1}->{ loess }=$tmpIQRMean;
    $loess->{-1}->{ stdev }=$tmpIQRStdev;
    
    if(validateLoessObject($loessObjectFile)) {
        open(COLLAPSED,outputWrapper($loessObjectFile,1)) or confess "Could not open file [$loessObjectFile] - $!";
        print COLLAPSED "-1\t-1\t$tmpIQRMean\t$tmpIQRStdev\n";
        print COLLAPSED "## done\n";
    } else {
        open(COLLAPSED,outputWrapper($loessObjectFile)) or confess "Could not open file [$loessObjectFile] - $!";
        print COLLAPSED "interactionDistance\trealInteractionDistance\tloessExpectedValue\tloessExpectedStdev\n";
        print COLLAPSED "-1\t-1\t$tmpIQRMean\t$tmpIQRStdev\n";
        print COLLAPSED "## done\n";
    }
    
    close(COLLAPSED);
    
    return($loess);
}

=head2 getRowColFactor

 Title     : getRowColFactor
 Usage     : $rowcolData=getRowColFactor(...)
 Function  : calculate row/col factor for normalization - dependent upon mode
 Returns   : primerData hash
 Argument  : matrixObject hash, input matrix file path, matrix name, includeCis flag, 
             maxDistance value, includeTrans flag, logTransform value, loess hash data hash, 
             skipZero flag, loess hash, factorMode, outputFile prefix, tmpFileName prefix
 
=cut

sub getRowColFactor($$$$$$$$$$$$;$$) {
    #required
    my $matrixObject=shift;
    my $inputMatrix=shift;
    my $inputMatrixName=shift;
    my $includeCis=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $includeTrans=shift;
    my $logTransform=shift;
    my $loess=shift;
    my $factorMode=shift;
    my $outputFile=shift;
    my $tmpFileName=shift;
    #optional 
    my $cisApproximateFactor=1;
    $cisApproximateFactor=shift if @_;
    my $excludeZero=0;
    $excludeZero=shift if @_;
    
    my $verbose=$matrixObject->{ verbose };
    
    confess "inputMatrix [$inputMatrix] does not exist." if(!(-e $inputMatrix));
    
    my $rowcolData={};

    my ($line);
    my ($lineNum,$dataIndex);
    $lineNum=$dataIndex=0;
    my %col2primer=();
   
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while($line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        $lineNum++;
        
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($dataIndex != 0) {
            my @tmp=split(/\t/,$line);
            my $anchorName=$tmp[0];
                        
            my $anchorObject={};
            $anchorObject=getHeaderObject($anchorName);
            
            my @tmpFactorArr=();
            my @tmpLoessMeanArr=();
            my @tmpLoessStdevArr=();
            
            my ($i);
            for($i=1;$i<@tmp;$i++) {
                my $partnerName=$col2primer{$i};
                my $partnerObject={};
                $partnerObject=getHeaderObject($partnerName);
                
                my $interactionDistance=getInteractionDistance($matrixObject,$anchorObject,$partnerObject,$cisApproximateFactor);
                my $interactionClassification=classifyInteraction($matrixObject,$includeCis,$minDistance,$maxDistance,$includeTrans,$anchorObject,$partnerObject);
                
                next if($interactionClassification ne "USABLE");
                
                my $cScore=$tmp[$i];
                
                next if($cScore eq "");
                next if(($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i));
                
                next if(($cScore == 0) and ($excludeZero));
                
                $cScore = (log($cScore+1)/log($logTransform)) if($logTransform != 0);
                
                confess "distance ($interactionDistance) does not exist!" if(!exists($loess->{$interactionDistance}));
                my $loessValue=$loess->{$interactionDistance}->{ loess };
                my $loessStdev=$loess->{$interactionDistance}->{ stdev };
                
                next if(($loessValue eq "NA") or ($loessStdev eq "NA")); 
                
                my $zScore = "NA";
                $zScore = (($cScore-$loessValue)/$loessStdev) if($loessStdev != 0);                
                my $obsExp = "NA";
                $obsExp = (log(($cScore/$loessValue))/log(2)) if(($loessValue ne "NA") and ($cScore > 0) and ($loessValue > 0));
                
                my $factor = "NA";
                $factor = $zScore if(($factorMode eq "zScore") or ($factorMode eq "zScore+obsExp"));
                $factor = $obsExp if($factorMode eq "obsExp");
                
                next if($factor eq "NA");
                
                push(@tmpFactorArr,$factor);
                push(@tmpLoessMeanArr,$loessValue);
                push(@tmpLoessStdevArr,$loessStdev);
                
            }
            
            my $tmpFactor = "NA";
            my $tmpLoessMean = "NA";
            my $tmpLoessStdev = "NA";
            
            if(@tmpFactorArr > 0) {
                my $tmpFactorArrStats=listStats(\@tmpFactorArr,0.10);
                $tmpFactor=$tmpFactorArrStats->{ iqrMean };
                my $tmpLoessMeanArrStats=listStats(\@tmpLoessMeanArr);
                $tmpLoessMean=$tmpLoessMeanArrStats->{ mean };
                my $tmpLoessStdevArrStats=listStats(\@tmpLoessStdevArr);
                $tmpLoessStdev=$tmpLoessStdevArrStats->{ mean };
            }
                        
            $tmpFactor = 2**$tmpFactor if(($factorMode eq "obsExp") and ($tmpFactor ne "NA"));
            $tmpFactor = sprintf "%.8f", $tmpFactor if($tmpFactor ne "NA");
            $rowcolData->{$anchorName}->{ factor }=$tmpFactor;    
            $rowcolData->{$anchorName}->{ loessMean }=$tmpLoessMean;    
            $rowcolData->{$anchorName}->{ loessStdev }=$tmpLoessStdev;    
            
        } else {
        
            #store the column headers index by column number
            my @tmp=split(/\t/,$line);
            
            my ($i);
            for($i=0;$i<@tmp;$i++) {
                $col2primer{$i}=$tmp[$i];
            }
            
        }
        
        $dataIndex++;
    }
    close(IN);
    
    open(OUT,outputWrapper($outputFile)) or confess "Could not open file [$outputFile] - $!";
    foreach my $primerName ( keys %$rowcolData ) {
        my $factor=$rowcolData->{$primerName}->{ factor };
        print OUT "$primerName\t$factor\n";
    }

    close(OUT);
    
    return($rowcolData);
}

=head2 getRowSum

 Title     : getRowSum
 Usage     : $rowSums=getRowSum(...)
 Function  : calculate sum of every row/col
 Returns   : rowSums hash
 Argument  : matrixObject hash, matrix 2D hash, maxDistance value, cisApproximateFactor value, 
             includeCis flag, includeTrans flag, ignoreZero flag, aggregrateMode value
 
=cut

sub getRowSum($$$$$$$$$) {
    # required
    my $matrixObject=shift;
    my $matrix=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $cisApproximateFactor=shift;
    my $includeCis=shift;
    my $includeTrans=shift;
    my $excludeZero=shift;
    my $aggregrateMode=shift;
    # optional 
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };

    my %rowSums=();
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        
        my @tmpList=();
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            
            #only work above diagonal if symmetrical 
            next if(($symmetrical) and ($y < $x)); 
            
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            
            my $cScore=$matrixObject->{ missingValue };
            $cScore=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
                        
            next if($cScore eq "");
            next if(($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i));
            next if($cScore eq "NA");
            
            next if(($excludeZero) and ($cScore == 0));
                        
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);
            my $interactionClassification=classifyInteraction($matrixObject,$includeCis,$minDistance,$maxDistance,$includeTrans,$yHeaderObject,$xHeaderObject);
                        
            next if($interactionClassification ne "USABLE");
            
            push(@tmpList,$cScore);
            
        }
        
        my $tmpArrStats=listStats(\@tmpList) if(@tmpList > 0);
        my $rowSum="NA";
        $rowSum=$tmpArrStats->{ $aggregrateMode } if(exists($tmpArrStats->{ $aggregrateMode }));
        
        $rowSums{$yHeader}=$rowSum;
    }
    
    return(\%rowSums);
}

=head2 matrix2listfile

 Title     : matrix2listfile
 Usage     : $sortedCisFile,$sortedTransFile,$matrixSum=matrix2listfile(...)
 Function  : dump matrix into cis/trans tsv files
 Returns   : cis file path, trans file path, matrix sum value
 Argument  : matrixObject hash, matrix 2D hash
 
=cut

sub matrix2listfile($$;$$$$$$$) {
    # required
    my $matrixObject=shift;
    my $inputMatrix=shift;
    # optional
    my $verbose=0;
    $verbose=shift if @_;
    my $excludeZero=0;
    $excludeZero=shift if @_;
    my $minDistance = undef;
    $minDistance = shift if @_;
    my $maxDistance = undef;
    $maxDistance=shift if @_;
    my $excludeCis=0;
    $excludeCis=shift if @_;
    my $includeCis=flipBool($excludeCis);
    my $excludeTrans=0;
    $excludeTrans=shift if @_;
    my $includeTrans=flipBool($excludeTrans);
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    print STDERR "\tmatrix2listfile\n\n" if($verbose);
        
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    my $headerSizing=$matrixObject->{ headerSizing };
    my $headerSpacing=$matrixObject->{ headerSpacing };
    my $missingValue=$matrixObject->{ missingValue };
    my $NA_rowcols=$matrixObject->{ NArowcols };
    my $output=$matrixObject->{ output };
    
    # create tmp dir
    my $tmpDir=createTmpDir();
    
    # cis file
    my $cisFile=$tmpDir.$inputMatrixName.".cis.txt.gz";
    open(CIS,outputWrapper($cisFile)) or confess "Could not open file [$cisFile] - $!" if($includeCis);
    
    # trans file
    my $transFile=$tmpDir.$inputMatrixName.".trans.txt.gz";
    open(TRANS,outputWrapper($transFile)) or confess "Could not open file [$transFile] - $!" if($includeTrans);
    
    my $lineNum=0;
    my @xHeaders=();
        
    my $nLines = getNumberOfLines($inputMatrix)-1;
    my $pcComplete=0;
    
    my $matrixSum=0;
    my $iCis=0;
    my $iTrans=0;
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            # x-headers
        } else {
            my @data=split(/\t/,$line);
            my $dsize=@data;
            my $yHeader=$data[0];
        
            # check if row is all NA
            next if(exists($NA_rowcols->{$yHeader}));
            
            my $yHeaderObject=getHeaderObject($yHeader);
            
            my $yIndex=-1;
            $yIndex = $header2inc->{ y }->{$yHeader} if(defined($header2inc->{ y }->{$yHeader}));
            next if($yIndex == -1);
            
            for(my $i=1;$i<@data;$i++) {
                my $cScore=$data[$i];
                
                # skip if cScore is not a valid number
                $cScore = "NA" if($cScore !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                
                next if($cScore eq "");
                next if(($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i));
                
                next if(($excludeZero) and ($cScore == 0));
                
                # truncate numbers to minimal digits
                $cScore = sprintf "%.".$sigDigits."f", $cScore if(($cScore ne "NA") and ($cScore !~ /^[+-]?\d+$/));
                
                my $xHeader="NA";
                $xHeader=$inc2header->{ x }->{$i-1};
                
                my $xIndex=-1;
                $xIndex = $header2inc->{ x }->{$xHeader} if(defined($header2inc->{ x }->{$xHeader}));
                next if($xIndex == -1);
                            
                my $xHeaderObject=getHeaderObject($xHeader);
                
                my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
                my $interactionClassification=classifyInteraction($matrixObject,$includeCis,$minDistance,$maxDistance,$includeTrans,$yHeaderObject,$xHeaderObject);

                next if($interactionClassification ne "USABLE");
                
                my $headerKey=$yHeader."___".$xHeader;
                
                if($interactionDistance == -1) {
                    next if($excludeTrans);
                    print TRANS "$interactionDistance\t$cScore\t$headerKey\n";
                    $matrixSum+=$cScore;
                    $iTrans++;
                } else { 
                    next if($excludeCis);
                    print CIS "$interactionDistance\t$cScore\t$headerKey\n";
                    $matrixSum+=$cScore;
                    $iCis++;
                }
            }
        }
        
        $pcComplete = 100 if($lineNum == ($nLines-1));
        print STDERR "\e[A" if($verbose);
        printf STDERR "\t%.2f%% complete ($lineNum/$nLines)...\n", $pcComplete if($verbose);
        $pcComplete = round((($lineNum/$nLines)*100),2);
        
        $lineNum++;
    }
    
    
    close(CIS) if($includeCis);
    close(TRANS) if($includeTrans);
    
    my $sortedCisFile=$output.".cis.txt.gz";
    system("gunzip -c '".$cisFile."' | sort -k1,1n | gzip >  '".$sortedCisFile."'") if($includeCis);
    system("rm '".$cisFile."'") if($includeCis);
    
    my $sortedTransFile=$output.".trans.txt.gz";
    system("gunzip -c '".$transFile."' | sort -k1,1n | gzip >  '".$sortedTransFile."'") if($includeTrans);
    system("rm '".$transFile."'") if($includeTrans);
    
    removeTmpDir($tmpDir);
    
    $pcComplete=100;
    print STDERR "\e[A" if($verbose);
    printf STDERR "\t%.2f%% complete ($lineNum/$nLines)...\n", $pcComplete if($verbose);
    
    return($sortedCisFile,$sortedTransFile,$matrixSum);
}

=head2 matrix2inputlist

 Title     : matrix2inputlist
 Usage     : $inputDataCis,$inputDataTrans=matrix2inputlist(...)
 Function  : matrix 2D hash to cis/trans arr ref
 Returns   : cis arr ref, trans arr ref
 Argument  : matrixObject hash, matrix 2D hash, includeCis flag, includeTrans flag
             maxDistance value, ignoreZero flag
 
=cut

sub matrix2inputlist($$$$$$$;$) {
    #required
    my $matrixObject=shift;
    my $matrix=shift;
    my $includeCis=shift;
    my $includeTrans=shift;
    my $minDistance=shift;
    my $maxDistance=shift;
    my $excludeZero=shift;
    #optional 
    my $cisApproximateFactor=1;
    $cisApproximateFactor=shift if @_;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    my $verbose=$matrixObject->{ verbose };
    
    my @inputDataTrans=();
    my @inputDataCis=();
    
    return(\@inputDataCis,\@inputDataTrans) if(($includeCis == 0) and ($includeTrans == 0));
    
    my $nInteractions=($numYHeaders*$numXHeaders)-1;
    my $progressBucketSize=1;
    $progressBucketSize=ceil($nInteractions/1000) if($nInteractions > 1000);
    my $pcComplete=0;
    
    print STDERR "\tmatrix2inputlist\n" if($verbose);
    
    my $iCis=0;
    my $iTrans=0;
    for(my $y=0;$y<$numYHeaders;$y++) {
        
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            
            #only work above diagonal if symmetrical 
            next if(($symmetrical) and ($y < $x)); 
            
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            
            my $cScore=$matrixObject->{ missingValue };
            $cScore=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            next if($cScore eq "");
            next if(($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i));
            
            next if(($excludeZero) and ($cScore == 0));
            
            my $headerKey=$yHeader."___".$xHeader;
            
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);
            
            # skip trans
            next if(($interactionDistance == -1) and ($includeTrans == 0));
            # skip cis
            next if(($interactionDistance != -1) and ($includeCis == 0));
            
            if($interactionDistance == -1) {
                $inputDataTrans[$iTrans][0]=$interactionDistance;
                $inputDataTrans[$iTrans][1]=$cScore;
                $inputDataTrans[$iTrans][2]=$headerKey;
                $iTrans++;
            } else { 
                my $realInteractionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1);
                $inputDataCis[$iCis][0]=$interactionDistance;
                $inputDataCis[$iCis][1]=$cScore;
                $inputDataCis[$iCis][2]=$yHeader."___".$xHeader;
                $inputDataCis[$iCis][3]=$realInteractionDistance;
                $iCis++;
            }
        }
        $pcComplete = round((($y/($numYHeaders-1))*100),2);
        print STDERR "\e[A" if(($verbose) and ($y != 0));
        printf STDERR "\t%.2f%% complete (".$y."/".($numYHeaders-1).")...\n", $pcComplete if($verbose);
    }
    
    # sort all input data by distance, smallest -> largest
    @inputDataCis = sort { $a->[3] <=> $b->[3] } @inputDataCis;        
    
    return(\@inputDataCis,\@inputDataTrans);
}

=head2 validateIdenticalMatrixStructure

 Title     : validateIdenticalMatrixStructure
 Usage     : $inc2header,$header2inc=validateIdenticalMatrixStructure(...)
 Function  : ensure two matrices have identitical structure (row/col headers)
 Returns   : inc2header hash, header2inc hash
 Argument  : matrix file 1 path, matrix file 2 path
             maxDistance value, ignoreZero flag
 
=cut

sub validateIdenticalMatrixStructure($$) {
    #required
    my $matrixFile_1=shift;
    my $matrixFile_2=shift;
    
    my ($inc2header_1,$header2inc_1)=parseHeaders($matrixFile_1);
    my $numYHeaders_1=keys(%{$inc2header_1->{ y }});
    my $numXHeaders_1=keys(%{$inc2header_1->{ x }});
    
    my ($inc2header_2,$header2inc_2)=parseHeaders($matrixFile_2);
    my $numYHeaders_2=keys(%{$inc2header_2->{ y }});
    my $numXHeaders_2=keys(%{$inc2header_2->{ x }});
    
    confess "#yFrags_1 ($numYHeaders_1) != #yFrags_2 ($numYHeaders_2)" if($numYHeaders_1 != $numYHeaders_2);
    confess "#xFrags_1 ($numXHeaders_1) != #xFrags_2 ($numXHeaders_2)" if($numXHeaders_1 != $numXHeaders_2);
    
    # check X Header overlaps
    for(my $x=0;$x<$numXHeaders_1;$x++) {
        my $xHeader_1=$inc2header_1->{ x }->{$x};
        $xHeader_1=deGroupHeader($xHeader_1,"liteChr",$x);
        
        my $xHeader_2="NA";
        $xHeader_2=$inc2header_2->{ x }->{$x} if(exists($inc2header_2->{ x }->{$x}));
        $xHeader_2=deGroupHeader($xHeader_2,"liteChr",$x);
        
        confess "xHeader_1 [$x] ($xHeader_1) != xHeader_2 [$x] ($xHeader_2)" if($xHeader_1 ne $xHeader_2);
    }
    for(my $x=0;$x<$numXHeaders_2;$x++) {
        my $xHeader_2=$inc2header_2->{ x }->{$x};
        $xHeader_2=deGroupHeader($xHeader_2,"liteChr",$x);
        
        my $xHeader_1="NA";
        $xHeader_1=$inc2header_1->{ x }->{$x} if(exists($inc2header_1->{ x }->{$x}));
        $xHeader_1=deGroupHeader($xHeader_1,"liteChr",$x);
        
        confess "xHeader_2 [$x] ($xHeader_2) != xHeader_1 [$x] ($xHeader_1)" if($xHeader_2 ne $xHeader_1);
    }
    
    # check Y Header overlaps
    for(my $y=0;$y<$numYHeaders_1;$y++) {
        my $yHeader_1=$inc2header_1->{ y }->{$y};
        $yHeader_1=deGroupHeader($yHeader_1,"liteChr",$y);
        
        my $yHeader_2="NA";
        $yHeader_2=$inc2header_2->{ y }->{$y} if(exists($inc2header_2->{ y }->{$y}));
        $yHeader_2=deGroupHeader($yHeader_2,"liteChr",$y);
        
        confess "yHeader_1 [$y] ($yHeader_1) != yHeader_2 [$y] ($yHeader_2)" if($yHeader_1 ne $yHeader_2);
    }
    for(my $y=0;$y<$numYHeaders_2;$y++) {
        my $yHeader_2=$inc2header_2->{ y }->{$y};
        $yHeader_2=deGroupHeader($yHeader_2,"liteChr",$y);
        
        my $yHeader_1="NA";
        $yHeader_1=$inc2header_1->{ y }->{$y} if(exists($inc2header_1->{ y }->{$y}));
        $yHeader_1=deGroupHeader($yHeader_1,"liteChr",$y);
        
        confess "yHeader_2 [$y] ($yHeader_2) != yHeader_1 [$y] ($yHeader_1)" if($yHeader_2 ne $yHeader_1);
    }
    
    #assume both inc2header_1 and inc2header_2 are identical, return one.
    
    return($inc2header_1,$header2inc_1);
    
}    

=head2 getHeaderSpacing

 Title     : getHeaderSpacing
 Usage     : $equalSpacingFlag,$equalSizingFlag,$meanGlobalHeaderSpacing,$meanGlobalHeaderSizing=getHeaderSpacing(...)
 Function  : calculatue bin spacing from array of headers
 Returns   : equalSpacing flag, equalSizing flag, averageHeaderSpacing, averageHeaderSize
 Argument  : inc2header sub-hash
 
=cut

sub getHeaderSpacing($) {
    my $inc2header=shift;
    
    my $numFrags=keys(%{$inc2header});
    
    my $equalSpacingFlag=1;
    my $equalSizingFlag=1;
    
    my (@globalHeaderSpacingArr,@globalHeaderSizingArr);
    my ($globalHeaderSpacing,$globalHeaderSizing);
    $globalHeaderSpacing=$globalHeaderSizing=-1;
    for(my $i=0;$i<$numFrags-1;$i++) {
        
        my $header=$inc2header->{$i};
        my $headerObject=getHeaderObject($header);
        my $headerRegion=$headerObject->{ region };
        my $headerStart=$headerObject->{ start };
        my $headerEnd=$headerObject->{ end };
        my $headerSize=$headerObject->{ size };
        
        my $nextHeader=$inc2header->{$i+1};
        my $nextHeaderObject=getHeaderObject($nextHeader);
        my $nextHeaderRegion=$nextHeaderObject->{ region };
        my $nextHeaderStart=$nextHeaderObject->{ start };
        my $nextHeaderEnd=$nextHeaderObject->{ end };
        my $nextHeaderSize=$nextHeaderObject->{ size };
        
        next if(($nextHeaderRegion ne $headerRegion) or ($headerEnd == $nextHeaderEnd) or ($headerStart == $nextHeaderStart));
        
        $equalSpacingFlag=0 if(($globalHeaderSpacing != (($nextHeaderStart-$headerStart))) and ($globalHeaderSpacing != -1));
        $equalSizingFlag=0 if(($globalHeaderSizing != ($headerSize)) and ($globalHeaderSizing != -1));
        
        $globalHeaderSpacing=($nextHeaderStart-$headerStart);
        $globalHeaderSizing=$headerSize;
                
        push(@globalHeaderSpacingArr,$globalHeaderSpacing);
        push(@globalHeaderSizingArr,$globalHeaderSizing);
        
    }
    
    my $meanGlobalHeaderSpacing=0;
    my $globalHeaderSpacingArrStats=listStats(\@globalHeaderSpacingArr) if(@globalHeaderSpacingArr > 0);
    $meanGlobalHeaderSpacing=$globalHeaderSpacingArrStats->{ mean } if(@globalHeaderSpacingArr > 0);
    
    
    my $meanGlobalHeaderSizing=0;
    my $globalHeaderSizingArrStats=listStats(\@globalHeaderSizingArr) if(@globalHeaderSizingArr > 0);
    $meanGlobalHeaderSizing=$globalHeaderSizingArrStats->{ mean } if(@globalHeaderSizingArr > 0);
    
    return($equalSpacingFlag,$equalSizingFlag,$meanGlobalHeaderSpacing,$meanGlobalHeaderSizing);
    
}

=head2 isSymmetrical

 Title     : isSymmetrical
 Usage     : $symmetrical=isSymmetrical(...)
 Function  : check if matrix headers are symmetrical
 Returns   : symmetrical
 Argument  : input matrix file OR inc2header hash
 
=cut

sub isSymmetrical($) {
    my $input=shift;
    
    # two possible inputs - either a file, or a hash ref to headers
    my ($inc2header,$header2inc);
    if(-e $input) {
        ($inc2header,$header2inc)=parseHeaders($input);
    } else {
        $inc2header=$input;
    }
        
    my $numYHeaders=keys(%{$inc2header->{ y }});
    my $numXHeaders=keys(%{$inc2header->{ x }});
    
    # lazy test for symmetrical heatmap
    return(0) if($numYHeaders != $numXHeaders); # not symmetrical

    # enforce perfectly symmetrical input matrix
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $xHeader=$inc2header->{ x }->{$y};
            
        return(0) if($yHeader ne $xHeader);
    }

    for(my $x=0;$x<$numXHeaders;$x++) {
        my $xHeader=$inc2header->{ x }->{$x};
        my $yHeader=$inc2header->{ y }->{$x};
        
        return(0) if($xHeader ne $yHeader);
    }
    
    # no longer do this 
    
    # if symmetrical headers, ensure header spacing is equal
    #my $numFrags=$numYHeaders;
    #my ($equalSpacingFlag_y,$equalSizingFlag_y,$headerSpacing_y,$headerSizing_y)=getHeaderSpacing($inc2header->{ y },$numYHeaders);
    #my ($equalSpacingFlag_x,$equalSizingFlag_x,$headerSpacing_x,$headerSizing_x)=getHeaderSpacing($inc2header->{ x },$numXHeaders);
    
    # enforce symmetrical headers, and all equal header spacing/sizing
    #return(0) if(($equalSpacingFlag_y == 0) or ($equalSizingFlag_y == 0)); # headers are not all equally sized/spaced
    #return(0) if(($equalSpacingFlag_x == 0) or ($equalSizingFlag_x == 0)); # headers are not all equally sized/spaced
    #return(0) if(($headerSpacing_y != $headerSpacing_x) or ($headerSizing_y != $headerSizing_x));
    
    return(1);
    
}

=head2 normalizeMatrix

 Title     : normalizeMatrix
 Usage     : $matrix=normalizeMatrix(...)
 Function  : normalize a matrix by read depth - scale by 10e6
 Returns   : matrix 2D hash
 Argument  : matrixObject hash, matrix 2D hash
 
=cut

sub normalizeMatrix($$;$$) {    
    my $matrixObject=shift;
    my $matrix=shift;
    #optional
    my $excludeDiagonal=0;
    $excludeDiagonal=shift if @_;
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    
    my $matrixSum=getMatrixSum($matrixObject,$matrix,$excludeDiagonal);
    
    my $sumMatrix=0;
    for(my $y=0;$y<$numYHeaders;$y++) {
        for(my $x=0;$x<$numXHeaders;$x++) {
        
            my $inten=$missingValue;
            $inten=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            next if($inten eq "NA");
            
            $inten = (($inten/$matrixSum)*1000000);
            $inten = sprintf "%.".$sigDigits."f", $inten;
            
            $matrix->{$y}->{$x}=$inten;
        }
    }
    
    return($matrix);
}

=head2 getMatrixSum

 Title     : getMatrixSum
 Usage     : $sum=getMatrixSum(...)
 Function  : return sum of matrix (non-symmetrical sum)
 Returns   : sum value
 Argument  : matrixObject hash, matrix 2D hash
 
=cut

sub getMatrixSum($$;$) {    
    my $matrixObject=shift;
    my $matrix=shift;
    #optional
    my $excludeDiagonal=0;
    $excludeDiagonal=shift if @_;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $missingValue=$matrixObject->{ missingValue };
    my $numFrags=$numYHeaders=$numXHeaders;
    
    my $symmetrical=isSymmetrical($inc2header);    
    
    my $sumMatrix=0;
    for(my $y=0;$y<$numYHeaders;$y++) {
        for(my $x=0;$x<$numXHeaders;$x++) {
            
            # skip below diagonal
            next if(($symmetrical) and ($y > $x));
            
            # skip diagonal
            next if(($excludeDiagonal) and ($y == $x) and ($symmetrical));
            
            my $inten=$missingValue;
            $inten=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
                        
            $sumMatrix += $inten if($inten ne "NA");
        }
    }
    
    return($sumMatrix);
}

=head2 getFileName

 Title     : getFileName
 Usage     : $fileName=getFileName(...)
 Function  : returns _stripped_ file name from file path
 Returns   : fileName
 Argument  : filePath string
 
=cut

sub getFileName($) {
    my $file=shift;
    
    my $fileName=(split(/\//,$file))[-1];
    my $shortName=$fileName;
    $shortName =~ s/\.matrix\.gz$//;
    $shortName =~ s/\.matrix$//;
    $shortName =~ s/\.gz$//;
    
    # if non-matrix file - remove extension
    $shortName=removeFileExtension($shortName) if($shortName eq $fileName);
    
    return($shortName);
}    

=head2 getShortFileName

 Title     : getShortFileName
 Usage     : $fileName=getShortFileName(...)
 Function  : returns _stripped_ and _shortened_ file name from file path
 Returns   : fileName
 Argument  : filePath string
 
=cut

sub getShortFileName($) {
    my $fileName=shift;
    
    $fileName=(split(/\//,$fileName))[-1];
    my $shortName=(split(/\./,$fileName))[0];
    $shortName=(split(/__/,$shortName))[0];
    
    return($shortName);
}    

=head2 getFilePath

 Title     : getFilePath
 Usage     : $filePath=getFilePath(...)
 Function  : returns filepath (name strimmed) - inverse of basename
 Returns   : filePath
 Argument  : filePath string
 
=cut

sub getFilePath($) {
    my $filePath=shift;
    
    my $shortName=(split(/\//,$filePath))[-1];
    $filePath =~ s/$shortName$//;    
    
    my $cwd = getcwd();
    $filePath = $cwd."/" if($filePath eq "");
    
    return($filePath);
}    

=head2 baseName

 Title     : baseName
 Usage     : $baseName=baseName(...)
 Function  : returns basename for file
 Returns   : baseName
 Argument  : filePath string
 
=cut

sub baseName($) {
    my $fileName=shift;
    
    my $shortName=(split(/\//,$fileName))[-1];
    
    return($shortName);
}    

=head2 removeFileExtension

 Title     : removeFileExtension
 Usage     : $baseName=removeFileExtension(...)
 Function  : strip file extension from file
 Returns   : fileName
 Argument  : filePath/fileName string
 
=cut

sub removeFileExtension($) {
    my $fileName=shift;
    
    my $extension=(split(/\./,$fileName))[-1];
    $fileName =~ s/\.$extension$//;
    
    return($fileName);
}

=head2 compareMatrices

 Title     : compareMatrices
 Usage     : $matrix=compareMatrices(...)
 Function  : compare two identical matrix 
 Returns   : matrix 2D hash
 Argument  : matrixObject_1 hash, matrixObject_2 hash, matrix_1 2D hash, matrix_2 2D hash, inc2header hash
 
=cut

sub compareMatrices($$$$$;$) {
    #required
    my $matrixObject_1=shift;
    my $matrixObject_2=shift;
    my $matrix_1=shift;
    my $matrix_2=shift;
    my $inc2header=shift;
    # optional
    my $compareMode="log2ratio";
    $compareMode=shift if @_;
    
    my $numYHeaders=keys(%{$inc2header->{ y }});
    my $numXHeaders=keys(%{$inc2header->{ x }});

    my %compareMatrix=();
    
    my ($y,$x);
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
        
            my $score_1=$matrixObject_1->{ missingValue };
            $score_1=$matrix_1->{$y}->{$x} if(defined($matrix_1->{$y}->{$x}));
            
            my $score_2=$matrixObject_2->{ missingValue };
            $score_2=$matrix_2->{$y}->{$x} if(defined($matrix_2->{$y}->{$x}));
            
            next if($score_1 eq "NA");
            next if($score_2 eq "NA");
            
            my $compareScore="NA";
            
            if($compareMode eq "log2ratio") {
                $compareScore=(log($score_1/$score_2)/log(2)) if(($score_1 > 0) and ($score_2 > 0));
            } elsif(($compareMode eq "add") or ($compareMode eq "sum")) {
                $compareScore=($score_1+$score_2);
            } elsif($compareMode eq "mean") {
                $compareScore=(($score_1+$score_2)/2);
            } elsif($compareMode eq "subtract") {
                $compareScore=($score_1-$score_2);
            } elsif($compareMode eq "divide") {
                $compareScore=($score_1/$score_2);
            } elsif($compareMode eq "multiply") {
                $compareScore=($score_1*$score_2);
            } elsif($compareMode eq "min") {
                $compareScore=min($score_1,$score_2);
            } elsif($compareMode eq "max") {
                $compareScore=max($score_1*$score_2);
            } else {
                confess "invalid compareMode [$compareMode]"
            }
        
            next if($compareScore eq "NA");
            $compareMatrix{$y}{$x}=$compareScore;
        }
    }
    
    return(\%compareMatrix);

}

=head2 correlateMatrices

 Title     : correlateMatrices
 Usage     : $matrix=correlateMatrices(...)
 Function  : correlate any two matrices, calculate scatter plot + R^2
 Returns   : tsv file, distance histogram hash
 Argument  : name, matrixObject_1 hash, matrixObject_2 hash, matrix_1 2D hash, matrix_2 2D hash, inc2header hash
 
=cut

sub correlateMatrices($$$$$$;$$$) {
    #required
    my $name=shift;
    my $matrixObject_1=shift;
    my $matrixObject_2=shift;
    my $matrix_1=shift;
    my $matrix_2=shift;
    my $inc2header=shift;
    # optional
    my $excludeZero=0;
    $excludeZero=shift if @_;
    my $excludeDiagonal=0;
    $excludeDiagonal=shift if @_;
    my $logTransform=0;
    $logTransform=shift if @_;
    
    my $correlationFile=$name.".correlate.txt.gz";
    open(OUT,outputWrapper($correlationFile)) or confess "Could not open file [$correlationFile] - $!";
    print OUT "interactionDistance\tcScore_1\tcScore_2\n";
    
    my $numYHeaders=keys(%{$inc2header->{ y }});
    my $numXHeaders=keys(%{$inc2header->{ x }});
    
    my $symmetrical_1=$matrixObject_1->{ symmetrical };
    my $symmetrical_2=$matrixObject_2->{ symmetrical };
    my $symmetrical = 0;
    $symmetrical = 1 if(($symmetrical_1) and ($symmetrical_2));
    
    my %distHist=();
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        for(my $x=0;$x<$numXHeaders;$x++) {
            # only work above diagonal if symmetrical map
            next if(($y > $x) and ($symmetrical)); 
            # skip diagonal if requested
            next if( ($excludeDiagonal) and (($symmetrical) and ($y == $x)) );
            
            my $score_1=$matrixObject_1->{ missingValue };
            $score_1=$matrix_1->{$y}->{$x} if(defined($matrix_1->{$y}->{$x}));
            next if($score_1 eq "NA");
            
            my $score_2=$matrixObject_2->{ missingValue };
            $score_2=$matrix_2->{$y}->{$x} if(defined($matrix_2->{$y}->{$x}));
            next if($score_2 eq "NA");
            
            next if( ($excludeZero) and (($score_1 == 0) and ($score_2 == 0)) );
            
            my $xHeader=$inc2header->{ x }->{$x};
            my $xHeaderObject=getHeaderObject($xHeader);
        
            my $interactionDistance=getInteractionDistance($matrixObject_1,$yHeaderObject,$xHeaderObject,1,$logTransform);
            
            print OUT "$interactionDistance\t$score_1\t$score_2\n";
            $distHist{$interactionDistance}++ if($interactionDistance > 0);
        }
    }
    
    close(OUT);
    
    return($correlationFile,\%distHist);

}

=head2 getMatrixAttributes

 Title     : getMatrixAttributes
 Usage     : $totalReads,$cisPercent,$transPercent,$averageTrans=getMatrixAttributes(...)
 Function  : return attributes for input matrix
 Returns   : totalReads [sum], cisPercent value, transPercent value, averageTrans value
 Argument  : input matrix file
 
=cut

sub getMatrixAttributes($;$) {
    # required
    my $inputMatrix=shift;
    # optional
    my $matrixObject={};
    $matrixObject=shift if @_;
    
    my $cisApproximateFactor=1;
    
    if(!exists($matrixObject->{ inc2header })) {
        $matrixObject=getMatrixObject($inputMatrix);
    }
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $symmetrical=$matrixObject->{ symmetrical };
    my $verbose=$matrixObject->{ verbose };
    
    my ($iCis,$iTrans,$totalCis,$totalTrans);
    $iCis=$iTrans=$totalCis=$totalTrans=0;
    
    my $lineNum=0;
    my @xHeaders=();
    
    print STDERR "\n\tgetMatrixAttributes\n" if($verbose);
    
    my $nLines = getNumberOfLines($inputMatrix);
    my $progressBucketSize=ceil($nLines / 1000);
    my $pcComplete=0;
    
    my $includeCis=1;
    my $includeTrans=1;
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            @xHeaders=split(/\t/,$line);
        } else {
            my @data=split(/\t/,$line);
            my $dsize=@data;
            my $yHeader=$data[0];
            
            my $yIndex=-1;
            $yIndex = $header2inc->{ y }->{$yHeader} if(defined($header2inc->{ y }->{$yHeader}));
            carp "header ($yHeader) does not exists in header2inc!" if($yIndex == -1);
            next if($yIndex == -1);
                
            my $yHeaderObject=getHeaderObject($yHeader);
            
            for(my $d=1;$d<$dsize;$d++) {
                my $xHeader=$xHeaders[$d];
                
                my $xIndex=-1;
                $xIndex = $header2inc->{ x }->{$xHeader} if(defined($header2inc->{ x }->{$xHeader}));
                carp "header ($xHeader) does not exists in header2inc!" if($xIndex == -1);
                next if($xIndex == -1);
            
                next if(($yIndex > $xIndex) and ($symmetrical)); # only work above diagonal if symmetrical map
    
                my $xHeaderObject=getHeaderObject($xHeader);            
                
                my $cScore=$data[$d];
                
                # skip if cScore is not a valid number
                next if($cScore !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                next if($cScore eq "");
                next if(($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i));
                
                my $interactionClassification=classifyInteraction($matrixObject,$includeCis,undef,undef,$includeTrans,$yHeaderObject,$xHeaderObject);
                next if($interactionClassification ne "USABLE");
                
                my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);
                
                if($interactionDistance == -1) {
                    $iTrans++;
                    $totalTrans += $cScore;
                } else { 
                    $iCis++;
                    $totalCis += $cScore;
                }
            }
        }
        
        $pcComplete = 100 if($lineNum == ($nLines-1));
        print STDERR "\e[A" if(($verbose) and ($lineNum != 0));
        printf STDERR "\t%.2f%% complete ($lineNum/$nLines)...\n", $pcComplete if($verbose);
        $pcComplete = round((($lineNum/$nLines)*100),2);
        $lineNum++;
        
    }
    close(IN);
    
    $pcComplete=100;
    print STDERR "\e[A" if($verbose);
    printf STDERR "\t%.2f%% complete ($lineNum/$nLines)...\n", $pcComplete if($verbose);
    
    my $totalReads=($totalCis+$totalTrans);
    
    my $cisPercent=round((($totalCis/$totalReads)*100),2);
    my $transPercent=round((($totalTrans/$totalReads)*100),2);
    
    my $averageTrans="NA";
    $averageTrans=($totalTrans/$iTrans) if($iTrans > 0);
    
    return($totalReads,$cisPercent,$transPercent,$averageTrans);
}

=head2 getMaxHeaderLength

 Title     : getMaxHeaderLength
 Usage     : $maxHeaderLength=getMaxHeaderLength(...)
 Function  : return maximum header length [in chars]
 Returns   : maxHeaderLength
 Argument  : inc2header hash
 
=cut

sub getMaxHeaderLength($) {
    my $inc2header=shift;
    
    my $numHeaders=keys(%{$inc2header});
    
    my $maxHeaderLength=0;
    for(my $i=0;$i<$numHeaders;$i++) {
        my $header=$inc2header->{$i};
        my $headerLength=length($header);
        $maxHeaderLength=$headerLength if($headerLength > $maxHeaderLength);
    }
    
    return($maxHeaderLength);
}

=head2 getMatrixObject

 Title     : getMatrixObject
 Usage     : $matrixObject=getMatrixObject(...)
 Function  : input matrix file to matrix object, collect all stats on matrix
 Returns   : matrixObject hash
 Argument  : input matrix file path
 
=cut

sub getMatrixObject($;$$$) {
    # required
    my $inputMatrix=shift;
    # optional
    my $tmpOutput="";
    $tmpOutput=shift if @_;
    my $verbose=0;
    $verbose=shift if @_;
    my $mode="";
    $mode=shift if @_;
    
    my $output=$tmpOutput;
    $output=getFileName($inputMatrix) if($output eq "");
    
    print STDERR "building matrix object [$mode]...\n" if($verbose);
    
    my %matrixObject=();
    
    # ensure file exists
    confess "file does not exist! [$inputMatrix]" if(!(-e $inputMatrix));
    
    # ensure file is valid
    confess "bad input file supplied (not normal matrix file) [$inputMatrix]" if(!(validateMatrixFile($inputMatrix)));
    
    # get matrix headers
    my ($headerFlag)=checkHeaders($inputMatrix);
    
    my ($inc2header,$header2inc)=parseHeaders($inputMatrix);
    my $numYHeaders=keys(%{$header2inc->{ y }});
    my $numXHeaders=keys(%{$header2inc->{ x }});
    my $numTotalHeaders=keys(%{$header2inc->{ xy }});
    
    # get matrix headers
    my ($header2contig,$index2contig,$contig2index,$contigList,$assembly)=parseContigs($inputMatrix,$inc2header,$header2inc);
    my $numXContigs=keys(%{$contig2index->{ x }});
    my $numYContigs=keys(%{$contig2index->{ y }});
    my $numContigs=keys(%{$contig2index->{ xy }});
    
    my $yMaxHeaderLength=getMaxHeaderLength($inc2header->{ y });
    my $xMaxHeaderLength=getMaxHeaderLength($inc2header->{ x });

    # get header spacing/sizing/overlap etc.
    my ($headerSizing,$headerSpacing,$binningStep,$equalSizingFlag,$equalSpacingFlag,$equalHeaderFlag)=getHeaderStats($inc2header);
    
    # check for matrix symmetry
    my $symmetrical=isSymmetrical($inc2header);
    
    # get NA row/cols
    my $NA_rowcols;
    if(($numTotalHeaders < 10000) and ($symmetrical) and ($mode ne "lite")) {
        my ($NA_rows)=getNARows($inputMatrix);
        my $NA_cols=$NA_rows;
        if($symmetrical == 0) {
            #print STDERR "transposing...\n";
            #my $transposedMatrix=transposeMatrix($inputMatrix);
            #($NA_cols)=getNARows($transposedMatrix);
            #system("rm '$transposedMatrix'");
        }
    
        my %NA_rowcols_hash=(%$NA_rows,%$NA_cols);
        $NA_rowcols=\%NA_rowcols_hash;
    }
    my $numNArowcols=keys %{$NA_rowcols};
    
    # calculate number of interactions
    my $numInteractions=($numYHeaders*$numXHeaders);
    confess "numInteractions [$numYHeaders x $numXHeaders] [$numInteractions] == 0!" if($numInteractions == 0);
    
    my $num_zeros="NA";
    my $num_nans="NA";
    my $num_nonzero_nonnan="NA";
    my $missingValue=0;
    ($num_zeros,$num_nans,$missingValue)=chooseOptimalMissingValue($inputMatrix) if(($numTotalHeaders < 20000) and ($mode ne "lite"));
    $num_nonzero_nonnan=($numInteractions-$num_zeros-$num_nans) if(($numTotalHeaders < 20000) and ($mode ne "lite"));
    
    my $inputMatrixName=getFileName($inputMatrix);
    my $inputMatrixPath=getFilePath($inputMatrix);
    
    if($verbose) {
        my $inputMatrixName_display=$inputMatrixName;
        $inputMatrixName_display=substr($inputMatrixName_display,0,60)."..." if(length($inputMatrixName_display) > 60);
        print STDERR "\tinputMatrixName\t$inputMatrixName_display\n";
        print STDERR "\tinputMatrixPath\t$inputMatrixPath\n";
        print STDERR "\tmatrixHeaderFlag\t$headerFlag\n";
        print STDERR "\tequalHeaderSizing\t$equalSizingFlag\n";
        print STDERR "\tequalHeaderSpacing\t$equalSpacingFlag\n";
        print STDERR "\theaderSizing\t$headerSizing\n";
        print STDERR "\theaderSpacing\t$headerSpacing\n";
        print STDERR "\tbinningStep\t$binningStep\n";
        print STDERR "\tequalHeaderFlag\t$equalHeaderFlag\n";
        print STDERR "\t# contigs\t".commify($numContigs)."\n";
        print STDERR "\t# yHeaders\t".commify($numYHeaders)."\n";
        print STDERR "\t# xHeaders\t".commify($numXHeaders)."\n";
        print STDERR "\t# totalHeaders\t".commify($numTotalHeaders)."\n";
        print STDERR "\tassembly\t$assembly\n";
        print STDERR "\tnumNArowcols\t".commify($numNArowcols)."\n";
        print STDERR "\tsymmetrical\t$symmetrical\n";
        print STDERR "\tnumInteractions\t".commify($numInteractions)."\n";
        print STDERR "\tclosest distance mode\n" if($equalHeaderFlag == 0);
        print STDERR "\tmidpoint distance mode\n" if($equalHeaderFlag);
        print STDERR "\tnum_zeros\t".commify($num_zeros)."\n";
        print STDERR "\tnum_nan\t".commify($num_nans)."\n";
        print STDERR "\tnum_nonzero_nonnan\t".commify($num_nonzero_nonnan)."\n";
        print STDERR "\tmissingValue\t$missingValue\n";
        print STDERR "\n";
    }
    
    $output=removeFileExtension($output) if($output =~ /\.gz$/);
    
    $matrixObject{ inputMatrixName }=$inputMatrixName;
    $matrixObject{ inputMatrixPath }=$inputMatrixPath;
    $matrixObject{ output }=$output;
    $matrixObject{ assembly }=$assembly;
    $matrixObject{ header2contig }=$header2contig;
    $matrixObject{ index2contig }=$index2contig;
    $matrixObject{ contig2index }=$contig2index;
    $matrixObject{ contigList }=$contigList;
    $matrixObject{ numXContigs }=$numXContigs;
    $matrixObject{ numYContigs }=$numYContigs;
    $matrixObject{ numContigs }=$numContigs;
    $matrixObject{ headerFlag }=$headerFlag;
    $matrixObject{ inc2header }=$inc2header;
    $matrixObject{ header2inc }=$header2inc;
    $matrixObject{ numYHeaders }=$numYHeaders;
    $matrixObject{ numXHeaders }=$numXHeaders;
    $matrixObject{ numTotalHeaders }=$numTotalHeaders;
    $matrixObject{ numNArowcols }=$numNArowcols;
    $matrixObject{ NArowcols }=$NA_rowcols;
    $matrixObject{ numInteractions }=$numInteractions;
    $matrixObject{ xHeaderLength }=$xMaxHeaderLength;
    $matrixObject{ yHeaderLength }=$yMaxHeaderLength;
    $matrixObject{ symmetrical }=$symmetrical;
    $matrixObject{ equalHeaderFlag }=$equalHeaderFlag;
    $matrixObject{ equalSizingFlag }=$equalSizingFlag;
    $matrixObject{ equalSpacingFlag }=$equalSpacingFlag;
    $matrixObject{ headerSizing }=$headerSizing;
    $matrixObject{ headerSpacing }=$headerSpacing;
    $matrixObject{ binningStep }=$binningStep;
    $matrixObject{ missingValue }=$missingValue;
    $matrixObject{ verbose }=$verbose;
    $matrixObject{ num_zeros }=$num_zeros;
    $matrixObject{ num_nans }=$num_nans;
    $matrixObject{ num_nonzero_nonnan }=$num_nonzero_nonnan;
    
    
    return(\%matrixObject);
    
}

=head2 getHeaderStats

 Title     : getHeaderStats
 Usage     : $headerSizing,$headerSpacing,$binningStep,$equalSizingFlag,$equalSpacingFlag,$equalHeaderFlag=getHeaderStats(...)
 Function  : collect information re. headers of input matrix
 Returns   : headerSizing, headerSpacing, binningStep, equalSizing flag, equalSpacing flag, equalHeader flagf
 Argument  : input matrix file path OR inc2header hash

=cut

sub getHeaderStats($) {
    my $input=shift;
    
    # two possible inputs - either a file, or a hash ref to headers
    my ($inc2header,$header2inc);
    if(-e $input) {
        ($inc2header,$header2inc)=parseHeaders($input);
    } else {
        $inc2header=$input;
    }
        
    my $numYHeaders=keys(%{$inc2header->{ y }});
    my $numXHeaders=keys(%{$inc2header->{ x }});
    
    my ($ySpacingFlag,$ySizingFlag,$yHeaderSpacing,$yHeaderSizing)=getHeaderSpacing($inc2header->{ y }); 
    my ($xSpacingFlag,$xSizingFlag,$xHeaderSpacing,$xHeaderSizing)=getHeaderSpacing($inc2header->{ x });

    my $equalSpacingFlag=0;
    $equalSpacingFlag = 1 if(($ySpacingFlag) and ($xSpacingFlag) and ($yHeaderSpacing == $xHeaderSpacing));
    
    my $equalSizingFlag=0;
    $equalSizingFlag = 1 if(($ySizingFlag) and ($xSizingFlag) and ($yHeaderSizing == $xHeaderSizing));
    
    my $equalHeaderFlag=0;
    my $binningStep=0;
    my $headerSpacing=-1;
    my $headerSizing=-1;
    if(($equalSpacingFlag) and ($equalSizingFlag)) {
        # y / x headers are equivalent
        $headerSpacing=$yHeaderSpacing=$xHeaderSpacing;
        $headerSizing=$yHeaderSizing=$xHeaderSizing;
        
        $binningStep=round($headerSizing/$headerSpacing) if($headerSpacing != 0);
        
        $equalHeaderFlag=1;
    }
    
    return($headerSizing,$headerSpacing,$binningStep,$equalSizingFlag,$equalSpacingFlag,$equalHeaderFlag);
    
}

=head2 correctMatrix

 Title     : correctMatrix
 Usage     : $matrix=correctMatrix(...)
 Function  : apply normalization to input matrix - derived from factors in primerData hash
 Returns   : matrix 2D hash
 Argument  : matrixObject hash, matrix 2D hash, primerData hash, includeCis flag, includeTrans flag, factorMode, logTransform, loess hash, ignoreZero flag

=cut

sub correctMatrix($$$$$$$$$;$) {
    #required
    my $matrixObject=shift;
    my $matrix=shift;
    my $rowcolData=shift;
    my $includeCis=shift;
    my $includeTrans=shift;
    my $factorMode=shift;
    my $logTransform=shift;
    my $loess=shift;
    my $excludeZero=shift;
    #optional
    my $cisApproximateFactor=1;
    $cisApproximateFactor=shift if @_;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $verbose=$matrixObject->{ verbose };
    
    my %normalizedMatrix=();
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader);
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};    
            my $xHeaderObject=getHeaderObject($xHeader);
            
            my $yFactor=$rowcolData->{$yHeader}->{ factor };
            my $xFactor=$rowcolData->{$xHeader}->{ factor };
            
            my $yLoessMean=$rowcolData->{$yHeader}->{ loessMean };
            my $xLoessMean=$rowcolData->{$xHeader}->{ loessMean };
            my $yLoessStdev=$rowcolData->{$yHeader}->{ loessStdev };
            my $xLoessStdev=$rowcolData->{$xHeader}->{ loessStdev };
            
            my $normalizedScore="NA";
            my $combinedFactor="NA";
    
            my $cScore=$matrixObject->{ missingValue };
            $cScore = $matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            $cScore = "NA" if(($cScore eq "") or ($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore eq "NA"));            
            
            # set to NA if cScore eq NA
            if($cScore eq "NA") {
                $normalizedMatrix{$y}{$x}="NA";
                next;
            }
            
            # set to NA if ignoreZero and cScore=0
            if(($cScore ne "NA") and (($excludeZero) and ($cScore == 0))) {
                $normalizedMatrix{$y}{$x}="NA";
                next;
            }
            
            if(($yFactor ne "NA") and ($xFactor ne "NA")) {
                
                my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,$cisApproximateFactor);    
                
                if( ($factorMode eq "zScore+obsExp") and ($cScore ne "NA")) {
                    
                    my $yCorrectionFactor = "NA";
                    $yCorrectionFactor = (($yLoessMean-($yLoessStdev*$yFactor))/($yLoessMean)) if(($yLoessMean ne "NA") and ($yLoessMean != 0));
                    my $xCorrectionFactor = "NA";
                    $xCorrectionFactor = (($xLoessMean-($xLoessStdev*$xFactor))/($xLoessMean)) if(($xLoessMean ne "NA") and ($xLoessMean != 0));
                    
                    if( (($yCorrectionFactor ne "NA") and ($yCorrectionFactor > 0)) and (($xCorrectionFactor ne "NA") and ($xCorrectionFactor > 0)) ) {
                        $combinedFactor = ($yCorrectionFactor*$xCorrectionFactor);
                        $normalizedScore = $cScore * $combinedFactor if($cScore ne "NA");
                    } else {
                        print STDERR "WARNING - ($yHeader | $xHeader) - $cScore [$yCorrectionFactor | $xCorrectionFactor] < 0!\n";    
                    }
                    
                } elsif( ($factorMode eq "zScore") and (($cScore ne "NA") and (exists($loess->{$interactionDistance}->{ loess }))) ) {
                    
                    $combinedFactor = ($yFactor+$xFactor);
                    my $expected = $loess->{$interactionDistance}->{ loess };
                    my $stdev = $loess->{$interactionDistance}->{ stdev };
                    $normalizedScore = ( $cScore - ( $combinedFactor*$loess->{$interactionDistance}->{ stdev } ) );
                    
                    $normalizedScore = 0 if(( $cScore - ( $combinedFactor*$loess->{$interactionDistance}->{ stdev } ) ) < 0);
                    #print STDERR "WARNING - ($yHeader | $xHeader\n\t$cScore [$yFactor | $xFactor | $combinedFactor] (".($cScore - ( $combinedFactor*$loess->{$interactionDistance}->{ stdev } ) )." < 0!)\n" if(( $cScore - ( $combinedFactor*$loess->{$interactionDistance}->{ stdev } ) ) <= 0);
                    
                } elsif($factorMode eq "obsExp") {
                    
                    $combinedFactor = ($yFactor*$xFactor);
                    $normalizedScore = $cScore * $combinedFactor if($cScore ne "NA");

                } else {    
                    print STDERR "ERROR! $factorMode | $cScore | $interactionDistance\n";
                }
                
            }
            
            $normalizedMatrix{$y}{$x}=$normalizedScore;
            
        }
    }
    
    return(\%normalizedMatrix);
}

=head2 validateZoomCoordinate

 Title     : validateZoomCoordinate
 Usage     : isValidZoom=validateZoomCoordinate(...)
 Function  : validate zoom coordinates (UCSC format)
 Returns   : flag
 Argument  : zoom coordinate string (UCSC format)

=cut

sub validateZoomCoordinate($) {
    my $zoomCoordinate=shift;
    
    $zoomCoordinate =~ s/,//g;
    
    return(0) if($zoomCoordinate eq "");
    return(0) if($zoomCoordinate !~ m/:/);
    return(0) if($zoomCoordinate !~ m/-/);
    
    my @tmp1=split(/:/,$zoomCoordinate);
    return(0) if(@tmp1 != 2);
    return(0) if($tmp1[0] eq "");
    
    my @tmp2=split(/-/,$tmp1[1]);
    return(0) if(@tmp2 != 2);
    return(0) if($tmp2[0] == $tmp2[1]);
    
    return(0) if($tmp2[0] !~ /^(\d+\.?\d*|\.\d+)$/);
    return(0) if($tmp2[1] !~ /^(\d+\.?\d*|\.\d+)$/);
    
    return(0) if($tmp2[0] > $tmp2[1]);
    
    return(1);
}

=head2 splitCoordinate

 Title     : splitCoordinate
 Usage     : coordinateData=splitCoordinate(...)
 Function  : split UCSC formatted coordinate into item hash
 Returns   : hash
 Argument  : coordinate string (UCSC format)

=cut

sub splitCoordinate($) {
    my $coordinate=shift;
    
    my ($coordinateChromosome,$coordinateStart,$coordinateEnd,$coordinatePosition);
    $coordinateChromosome=$coordinateStart=$coordinateEnd=$coordinatePosition="NA";
    
    my $goodCoordinateFlag=0;
    if(($coordinate ne "") and (validateZoomCoordinate($coordinate))) {
    
        ($coordinateChromosome,$coordinatePosition)=split(/:/,$coordinate);
        $coordinatePosition =~ s/,//g;
    
        ($coordinateStart,$coordinateEnd)=split(/-/,$coordinatePosition);
        
        $goodCoordinateFlag=1;
    }
     
    my %coordinateData=();
    $coordinateData{ chromosome } = $coordinateChromosome;
    $coordinateData{ start } = $coordinateStart;
    $coordinateData{ end } = $coordinateEnd;
    $coordinateData{ size } = ($coordinateEnd-$coordinateStart);
    $coordinateData{ flag } = $goodCoordinateFlag;
    $coordinateData{ name } = $coordinateChromosome.":".$coordinateStart."-".$coordinateEnd;
    
    return(\%coordinateData);
}

=head2 header2subMatrix

 Title     : header2subMatrix
 Usage     : $subMatrix=header2subMatrix(...)
 Function  : return sub-group of header [region, chr, nakedChr, liteChr, group]
 Returns   : string
 Argument  : header string, extractBy string

=cut

sub header2subMatrix($$) {
    my $header=shift;
    my $extractBy=shift;
    
    my $headerObject=getHeaderObject($header);
    
    my $chromosome=$headerObject->{ chromosome };
    my $region=$headerObject->{ region };
        
    my @tmp=split(/-/,$chromosome);
    my $group="amb";
    $group=$tmp[1] if(@tmp == 2);
    
    my $nakedChromosome=$chromosome;
    $nakedChromosome =~ s/chr//;
    
    my $liteChromosome=$chromosome;
    $liteChromosome =~ s/-$group//;
    
    my $subMatrix="NA";
    $subMatrix=$region if($extractBy eq "region");
    $subMatrix=$chromosome if($extractBy eq "chr");
    $subMatrix=$nakedChromosome if($extractBy eq "nakedChr");
    $subMatrix=$liteChromosome if($extractBy eq "liteChr");
    $subMatrix=$group if($extractBy eq "group");
    
    return($subMatrix);
}

=head2 deGroupHeader

 Title     : deGroupHeader
 Usage     : $header=deGroupHeader(...)
 Function  : strip off allelic info from chromosome
 Returns   : string
 Argument  : header string

=cut

sub deGroupHeader($;$$) {
    #required 
    my $header=shift;
    #optional
    my $extractBy="liteChr";
    $extractBy=shift if @_;
    my $index=undef;
    $index=shift if @_;
    $index=getSmallUniqueString() if(!defined($index));
    
    my $headerObject=getHeaderObject($header);

    my $subName=$headerObject->{ subName };
    my $assembly=$headerObject->{ assembly };
    my $chromosome=$headerObject->{ chromosome };
    my $start=$headerObject->{ start };
    my $end=$headerObject->{ end };

    my $liteChromosome=header2subMatrix($header,$extractBy);

    my $deGroupedHeader=$index."|".$assembly."|".$liteChromosome.":".$start."-".$end;
    
    return($deGroupedHeader);
}

=head2 reOrientIntervals

 Title     : reOrientIntervals
 Usage     : $start1,$end1,$start2,$end2=reOrientIntervals(...)
 Function  : re-orient two intervals, so that interval_a < interval_b
 Returns   : start_pos_a, end_pos_a, start_pos_b, end_pos_b
 Argument  : start_pos_a, end_pos_a, start_pos_b, end_pos_b

=cut

sub reOrientIntervals($$$$) {
    my $start1=shift;
    my $end1=shift;
    my $start2=shift;
    my $end2=shift;
    
    return($start1,$end1,$start2,$end2) if($start1 <= $start2);
    return($start2,$end2,$start1,$end1) if($start2 <= $start1);
    
    confess "poorly formed intervals! $start1 - $end1 | $start2 - $end2";
}

=head2 isOverlapping

 Title     : isOverlapping
 Usage     : $start1,$end1,$start2,$end2=isOverlapping(...)
 Function  : check if two intervals overlap
 Returns   : flag
 Argument  : start_pos_a, end_pos_a, start_pos_b, end_pos_b

=cut

sub isOverlapping($$$$;$$) {
    #required
    my $start1=shift;
    my $end1=shift;
    my $start2=shift;
    my $end2=shift;
    #optional
    my $chromosome1="NA";
    $chromosome1=shift if @_;
    my $chromosome2="NA";
    $chromosome2=shift if @_;
    
    return(0) if($chromosome1 ne $chromosome2);
    
    confess "poorly formed intervals! $start1 - $end1 | $start2 - $end2" if(($start1 > $end1) or (!defined($start1)) or (!defined($end1)));
    confess "poorly formed intervals! $start1 - $end1 | $start2 - $end2" if(($start2 > $end2) or (!defined($end1)) or (!defined($end2)));
    
    ($start1,$end1,$start2,$end2)=reOrientIntervals($start1,$end1,$start2,$end2);

    # method 1 
    my $overlap_1 = 0;
    $overlap_1 = 1 if(max($start1,$start2) <= min($end1,$end2));

    # method 2
    my $overlap_2 = 0;
    $overlap_2 = 1 if(($start1 <= $end2) and ($start2 <= $end1));
        
    return(1) if(($overlap_1) and ($overlap_2));
    
    confess "overlap logic failire ($overlap_1 || $overlap_2)" if($overlap_1 != $overlap_2);
    
    return(0);
}

=head2 stripChromosomeGroup

 Title     : stripChromosomeGroup
 Usage     : $chr=stripChromosomeGroup(...)
 Function  : remove allelic group from chromosome name 
 Returns   : string
 Argument  : chromosome string
 Example   : chr1-129S1 -> chr1

=cut

sub stripChromosomeGroup($) {
    my $chromosome=shift;
    
    my @tmp=split(/-/,$chromosome);
    my $group="amb";
    $group=$tmp[1] if(@tmp == 2);
    
    my $liteChromosome=$chromosome;
    $liteChromosome =~ s/-$group//;
    
    return($liteChromosome);
}

=head2 stripAssembly

 Title     : stripAssembly
 Usage     : $assembly=stripAssembly(...)
 Function  : remove allelic group from assembly name 
 Returns   : string
 Argument  : assembly
 Example   : mm9-cast-129s1 -> mm9

=cut

sub stripAssemblyGroup($) {
    my $assembly=shift;
    
    my @tmp=split(/-/,$assembly);
    my $liteAssembly=$assembly;
    $liteAssembly=$tmp[0] if(@tmp == 3);

    return($liteAssembly);
}

=head2 getRestrictionEnzymeSequences

 Title     : getRestrictionEnzymeSequences
 Usage     : $rs=getRestrictionEnzymeSequences()
 Function  : build hash of enzyme/cut-sequences
 Returns   : hash
 Argument  : 

=cut

sub getRestrictionEnzymeSequences() {
    my %restrictionEnzymeSequences=();
    
    $restrictionEnzymeSequences{ HindIII } = "AAGCTT";
    $restrictionEnzymeSequences{ EcoRI } = "GAATTC";
    $restrictionEnzymeSequences{ NcoI } = "CCATGG";
    $restrictionEnzymeSequences{ DpnII } = "GATC";
    $restrictionEnzymeSequences{ MNase } = "MNase";
    
    return(\%restrictionEnzymeSequences);
}
    
=head2 getUserHomeDirectory

 Title     : getUserHomeDirectory
 Usage     : $home=getUserHomeDirectory()
 Function  : return user home dir
 Returns   : file path
 Argument  : 

=cut

sub getUserHomeDirectory() {
    my $userHomeDirectory = `echo \$HOME`;
    chomp($userHomeDirectory);
    return($userHomeDirectory);
}

=head2 getUniqueString

 Title     : getUniqueString
 Usage     : $uniq=getUniqueString()
 Function  : return randomly generated unique string [via UUID]
 Returns   : string
 Argument  : 

=cut

sub getUniqueString() {
    my $UUID = `uuidgen`;
    chomp($UUID);
    return($UUID);
}

=head2 getSmallUniqueString

 Title     : getSmallUniqueString
 Usage     : $uniq=getSmallUniqueString()
 Function  : return randomly generated _small_ unique string [via UUID]
 Returns   : string
 Argument  : 

=cut

sub getSmallUniqueString() {
    my $UUID=`uuidgen | rev | cut -d '-' -f 1`;
    chomp($UUID);
    return($UUID);
}

=head2 getComputeResource

 Title     : getComputeResource
 Usage     : $host=getComputeResource()
 Function  : return host name 
 Returns   : string
 Argument  : 

=cut

sub getComputeResource() {
    my $hostname = `hostname`;
    chomp($hostname);
    return($hostname);
}

=head2 translateFlag

 Title     : translateFlag
 Usage     : $flag=translateFlag(...)
 Function  : translate bool flag to on/off verbage
 Returns   : string
 Argument  : flag bool

=cut

sub translateFlag($) {
    my $flag=shift;
    
    my $response="off";
    $response="on" if($flag);    
    return($response);
}

=head2 getNumberOfLines

 Title     : getNumberOfLines
 Usage     : $nLines=getNumberOfLines(...)
 Function  : get number of lines in file
 Returns   : number
 Argument  : input file path

=cut

sub getNumberOfLines($) {
    my $inputFile=shift;
    
    my $nLines = 0;
    
    if(($inputFile =~ /\.gz$/) and (!(-T($inputFile)))) {
        my $matrixInfo=`gunzip -c '$inputFile' 2>/dev/null | head -n 1 | cut -f 1`;
        if($matrixInfo =~ /(\d+)x(\d+)/) {
            chomp($matrixInfo);
            my ($nCols,$nRows) = split(/x/,$matrixInfo);
            $nLines = $nRows;
        } else { 
            $nLines = `gunzip -c '$inputFile' 2>/dev/null | grep -v "^#" | grep -v "^\$" | wc -l`;
        }
    } else {
        $nLines = `grep -v "^#" '$inputFile' | grep -v "^\$" | wc -l`;
    }

    chomp($nLines);
    $nLines =~ s/ //g;
    
    return($nLines);
}
    
=head2 matrix2pairwise

 Title     : matrix2pairwise
 Usage     : $pairfile=matrix2pairwise(...)
 Function  : transform matrix into 3 colulmn tsv file (header1, header2, score)
 Returns   : file path
 Argument  : matrixObject hash, input matrix file path, matrix name

=cut

sub matrix2pairwise($$;$$$$$) {
    # required
    my $matrixObject=shift;
    my $inputMatrix=shift;
    # optional
    my $excludeCis=0;
    $excludeCis=shift if @_;
    my $excludeTrans=0;
    $excludeTrans=shift if @_;
    my $skipNAs=0;
    $skipNAs=shift if @_;
    my $excludeZero=0;
    $excludeZero=shift if @_;
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    my $verbose=$matrixObject->{ verbose };
    my $output=$matrixObject->{ output };
    
    my $tmpDir=createTmpDir();
    my $tmpPairwiseFile=$tmpDir.$output.".pairwise.txt.gz";
    open(OUT,outputWrapper($tmpPairwiseFile)) or confess "Could not open file [$tmpPairwiseFile] - $!";
    
    my $lineNum=0;
    my @xHeaders=();
    
    my $nLines = getNumberOfLines($inputMatrix)-2;
    my $progressBucketSize=ceil($nLines / 1000);
    my $pcComplete=0;
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            @xHeaders=split(/\t/,$line);
        } else {
            my @data=split(/\t/,$line);
            my $dsize=@data;
            my $yHeader=$data[0];
            my $yHeaderObject=getHeaderObject($yHeader);
            
            for(my $d=1;$d<$dsize;$d++) {
                my $xHeader=$xHeaders[$d];
                my $xHeaderObject=getHeaderObject($xHeader);

                my $cScore=$data[$d];
                
                $cScore = "NA" if($cScore !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                $cScore = "NA" if(($cScore eq "") or ($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i) or ($cScore eq "NA"));
                
                next if(($skipNAs) and ($cScore eq "NA")); # skip NAs if selected
                next if($cScore ne "NA") and (($excludeZero) and ($cScore == 0));
                
                $cScore = sprintf "%.".$sigDigits."f", $cScore if($cScore ne "NA");
                
                my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject);
                
                next if(($interactionDistance == -1) and ($excludeTrans));
                next if(($interactionDistance != -1) and ($excludeCis));
                
                print OUT "$yHeader\t$xHeader\t$cScore\n";
                
            }
        }
        $pcComplete = 100 if($lineNum >= $nLines);
        print STDERR "\e[A" if(($verbose) and ($lineNum != 0));
        printf STDERR "\t%.2f%% complete ($lineNum/".$nLines.")...\n", $pcComplete if($verbose);
        $pcComplete = round((($lineNum/$nLines)*100),2);
        $lineNum++;
    }
    close(IN);
    
    close(OUT);
    
    my $tmpCommentFile=$tmpDir.$inputMatrixName.".comments.txt";
    
    # get orig comment lines
    system("gunzip -c '".$tmpPairwiseFile."' | grep -v '^[^#;]' > '".$tmpCommentFile."'");
    
    # append on comment headers
    open(COMMENT,outputWrapper($tmpCommentFile,0,1)) or confess "Could not open file [$tmpCommentFile] - $!";
    print COMMENT "yHeader\txHeader\tcScore\n";
    close(COMMENT);
    
    # gzip comment file for cat
    system("gzip '".$tmpCommentFile."'");
    
    # do the sort, remove comment lines
    my $sortedPairwiseFile=$tmpDir.$inputMatrixName.".sorted.pairwise-score.txt.gz";
    system("gunzip -c '".$tmpPairwiseFile."' | grep '^[^#;]' | sort -k3,3n | gzip > '".$sortedPairwiseFile."'");
    
    # put comment lines back on
    my $pairwiseFile=$output.".pairwise-score.txt.gz";
    system("cat '".$tmpCommentFile.".gz' '".$sortedPairwiseFile."' > '".$pairwiseFile."'");
    
    removeTmpDir($tmpDir);
    
    return($pairwiseFile);
}

=head2 matrix2distance

 Title     : matrix2distance
 Usage     : $distfile=matrix2distance(...)
 Function  : transform matrix into 3 colulmn tsv file (header1, header2, distance)
 Returns   : file path
 Argument  : matrixObject hash, input matrix file path, matrix name

=cut

sub matrix2distance($$;$$$$$) {
    # required
    my $matrixObject=shift;
    my $inputMatrix=shift;
    # optional
    my $excludeCis=0;
    $excludeCis=shift if @_;
    my $excludeTrans=0;
    $excludeTrans=shift if @_;
    my $skipNAs=0;
    $skipNAs=shift if @_;
    my $excludeZero=0;
    $excludeZero=shift if @_;
    my $sigDigits=4;
    $sigDigits=shift if @_;
    
    my $verbose=$matrixObject->{ verbose };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    my $output=$matrixObject->{ output };
    
    my $tmpDir=createTmpDir();
    my $tmpPairwiseFile=$tmpDir.$inputMatrixName.".pairwise-distance.txt.gz";
    open(OUT,outputWrapper($tmpPairwiseFile)) or confess "Could not open file [$tmpPairwiseFile] - $!";
    
    my $lineNum=0;
    my @xHeaders=();
    
    my $nLines = getNumberOfLines($inputMatrix)-2;
    my $progressBucketSize=ceil($nLines / 1000);
    my $pcComplete=0;
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum == 0) {
            @xHeaders=split(/\t/,$line);
        } else {
            my @data=split(/\t/,$line);
            my $dsize=@data;
            my $yHeader=$data[0];
            my $yHeaderObject=getHeaderObject($yHeader);
            
            for(my $d=1;$d<$dsize;$d++) {
                my $xHeader=$xHeaders[$d];
                my $xHeaderObject=getHeaderObject($xHeader);

                my $cScore=$data[$d];
                $cScore = "NA" if($cScore !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
                $cScore = "NA" if(($cScore eq "") or ($cScore =~ /^NULL$/i) or ($cScore =~ /^NA$/i) or ($cScore =~ /inf$/i) or ($cScore =~ /^nan$/i) or ($cScore eq "NA"));
                
                next if(($skipNAs) and ($cScore eq "NA")); # skip NAs if selected
                next if($cScore ne "NA") and (($excludeZero) and ($cScore == 0));
                $cScore = sprintf "%.".$sigDigits."f", $cScore if($cScore ne "NA");
                
                my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject);
                
                next if(($interactionDistance == -1) and ($excludeTrans));
                next if(($interactionDistance != -1) and ($excludeCis));
                
                print OUT "$yHeader\t$xHeader\t$interactionDistance\n";
                
            }
        }
        $pcComplete = 100 if($lineNum >= $nLines);
        print STDERR "\e[A" if(($verbose) and ($lineNum != 0));
        printf STDERR "\t%.2f%% complete ($lineNum/".$nLines.")...\n", $pcComplete if($verbose);
        $pcComplete = round((($lineNum/$nLines)*100),2);
        $lineNum++;
    }
    close(IN);
    
    close(OUT);
    
    
    my $tmpCommentFile=$tmpDir.$inputMatrixName.".comments.txt";
    
    # get orig comment lines
    system("gunzip -c '".$tmpPairwiseFile."' | grep -v '^[^#;]' > '".$tmpCommentFile."'");
    
    # append on comment headers
    open(COMMENT,outputWrapper($tmpCommentFile,0,1)) or confess "Could not open file [$tmpCommentFile] - $!";
    print COMMENT "yHeader\txHeader\tdistance\n";
    close(COMMENT);
    
    # gzip comment file for cat
    system("gzip '".$tmpCommentFile."'");
    
    # do the sort, remove comment lines
    my $sortedPairwiseFile=$tmpDir.$inputMatrixName.".sorted.pairwise-distance.txt.gz";
    system("gunzip -c '".$tmpPairwiseFile."' | grep '^[^#;]' | sort -k3,3n | gzip > '".$sortedPairwiseFile."'");
    
    # put comment lines back on
    my $pairwiseFile=$output.".pairwise-distance.txt.gz";
    system("cat '".$tmpCommentFile.".gz' '".$sortedPairwiseFile."' > '".$pairwiseFile."'");
    
    removeTmpDir($tmpDir);
    
    return($pairwiseFile);
}

=head2 autoSize

 Title     : autoSize
 Usage     : $pixelSize=autoSize(...)
 Function  : calculate nice pixelSize for image given width/height
 Returns   : pixel size
 Argument  : ideal heatmap size [pixels], num y-axis headers, num x-axis headers

=cut

sub autoSize($$$) {
    my $heatmapSize=shift;
    my $numYHeaders=shift;
    my $numXHeaders=shift;
    
    my ($pixelSize,$imageHeight,$imageWidth);
    
    if($numXHeaders > $numYHeaders) {
        $pixelSize=ceil($heatmapSize/$numXHeaders);
    } else {
        $pixelSize=ceil($heatmapSize/$numYHeaders);
    }
    
    return($pixelSize);
}

=head2 autoScale

 Title     : autoScale
 Usage     : $start,$end=autoScale(...)
 Function  : calculate nice color-range for array list of signals
 Returns   : start,end
 Argument  : matrixObject hash, matrix 2D hash

=cut

sub autoScale($$;$$$) {
    #required
    my $matrixObject=shift;
    my $matrix=shift;
    # optional
    my $startQuantile=0.025;
    $startQuantile=shift if @_;
    my $endQuantile=0.975;
    $endQuantile=shift if @_;
    my $interactionType="all";
    $interactionType=shift if @_;
    
    $interactionType="all" if(($interactionType ne "cis") and ($interactionType ne "trans"));
    
    my $inc2header=$matrixObject->{ inc2header };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $symmetrical = $matrixObject->{ symmetrical };
    my $missingValue = $matrixObject->{ missingValue };
    my $num_zeros = $matrixObject->{ num_zeros };
    my $num_interactions=$matrixObject->{ numInteractions };    
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    my $output=$matrixObject->{ output };
    
    my $tmpDir=createTmpDir();
    my $tmpFile=$tmpDir.$inputMatrixName.".autoscale.gz";
    open(OUT,outputWrapper($tmpFile)) or confess "Could not open file [$tmpFile] - $!";
    
    my $widthPC=0;
    $widthPC=0.02 if($symmetrical);
    my $diagonalWidth=ceil($numYHeaders*$widthPC);
    # only subset exact diagonal 
    $diagonalWidth=1;
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        my $yHeader=$inc2header->{ y }->{$y};
        my $yHeaderObject=getHeaderObject($yHeader) if($interactionType ne "all");
        
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $xHeader=$inc2header->{ x }->{$x};
            
            # skip below diagonal
            next if(($symmetrical) and ($y > $x));
            
            next if($yHeader eq $xHeader);
            next if(($symmetrical) and (abs($y-$x) <= $diagonalWidth));
            
            my $score = $missingValue;
            $score=$matrix->{$y}->{$x} if(defined($matrix->{$y}->{$x}));
            
            next if($score eq "NA");
            next if( ($score eq $missingValue) or (($score ne "NA") and ($missingValue ne "NA") and ($score == $missingValue)) );
            
            my $xHeaderObject=getHeaderObject($xHeader) if($interactionType ne "all");
            my $interactionDistance=getInteractionDistance($matrixObject,$yHeaderObject,$xHeaderObject,1) if($interactionType ne "all");
            
            if($interactionType ne "all") {
                next if(($interactionType eq "cis") and ($interactionDistance == -1));
                next if(($interactionType eq "trans") and ($interactionDistance != -1));
            }
            
            $score=abs($score);
            
            print OUT "$score\n"
            
        }
    }
    
    close(OUT);
    
    my $colorScaleStart=0;
    my $colorScaleEnd=1;
    
    my $sortedTmpFile=$tmpDir.$inputMatrixName.".autoscale.sorted.gz";
    system("gunzip -c '".$tmpFile."' | sort -k1,1n | gzip >  '".$sortedTmpFile."'");
    
    my $nLines = getNumberOfLines($sortedTmpFile);
    
    my $bottomLine = ceil($startQuantile * $nLines);
    my $topLine = floor($endQuantile * $nLines);
    
    my $lineNum=0;
    
    open(IN,inputWrapper($sortedTmpFile)) or confess "Could not open file [$sortedTmpFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        $lineNum++;
        next if($lineNum < $bottomLine); # skip beginning
        next if(($lineNum > $bottomLine) and ($lineNum < $topLine)); # skip middle
        last if($lineNum > $topLine); # skip end
        
        $colorScaleStart=$line if($lineNum == $bottomLine);
        $colorScaleEnd=$line if($lineNum == $topLine);
    }
    close(IN);

    system("rm '".$tmpFile."'");
    system("rm '".$sortedTmpFile."'");
    
    $colorScaleStart=0 if(($missingValue eq 0) and (($num_zeros/$num_interactions) > $startQuantile));
    
    if($colorScaleEnd <= 0) { $colorScaleEnd = 1; }
    my $orig_colorScaleStart=$colorScaleStart;
    if($colorScaleStart == $colorScaleEnd) { 
        $colorScaleStart=($colorScaleEnd-1); 
        $colorScaleEnd=$orig_colorScaleStart+1;
    }
    
    removeTmpDir($tmpDir);

    return($colorScaleStart,$colorScaleEnd);
}    

=head2 classifyInteractionDistance

 Title     : classifyInteractionDistance
 Usage     : $classification=classifyInteractionDistance(...)
 Function  : classify interaction by distance
 Returns   : classification
 Argument  : distance

=cut

sub classifyInteractionDistance($) {
    my $interactionDistance=shift;
    
    return("cis") if($interactionDistance != -1);
    return("trans") if($interactionDistance == -1);
    return("NA")
}

=head2 getAvailableColorScales

 Title     : getAvailableColorScales
 Usage     : $availableColors,$transparency=getAvailableColorScales()
 Function  : generate all color scales, dump to hash
 Returns   : availableColor hash, transparency value
 Argument  : 

=cut

sub getAvailableColorScales() {
    
    my %availableColorScales=();
    
    # rainbow color scale
    $availableColorScales{ rainbow } = "red,orange,yellow,green,blue,indigo,violet"; 
    # inverse rainbow color scale
    $availableColorScales{ inverseRainbow } = "violet,indigo,blue,green,yellow,orange,red"; 
    
    # dekker standard color scales
    $availableColorScales{ dekkerPositive } = "white,orange,red,darkRed"; 
    $availableColorScales{ dekkerNegative } = "white,cyan,blue,darkBlue"; 
    
    # dekker standard color scales
    $availableColorScales{ americaPositive } = "white,lightRed,red,darkRed"; 
    $availableColorScales{ americaNegative } = "white,lightBlue,blue,darkBlue"; 
    
    
    return(\%availableColorScales);
}

=head2 getAvailableColors

 Title     : getAvailableColors
 Usage     : $availableColors,$transparency=getAvailableColors()
 Function  : generate all available RGB colors
 Returns   : availableColor hash, transparency value
 Argument  : 

=cut

sub getAvailableColors() {
    
    my $transparency=0;
    my %availableColors=();
    
    # image background
    @{$availableColors{ white }} = (255,255,255,$transparency); 
    
    # utility colorPalette
    @{$availableColors{ black }} = (0,0,0,$transparency);
    @{$availableColors{ null }} = (150,150,150,0);
    @{$availableColors{ grey }} = (100,100,100,$transparency);
    #primary
    @{$availableColors{ red }} = (255,0,0,$transparency);
    @{$availableColors{ green }} = (0,255,0,$transparency);
    @{$availableColors{ blue }} = (0,0,255,$transparency);
    #secondary
    @{$availableColors{ yellow }} = (255,255,0,$transparency);
    @{$availableColors{ violet }} = (255,0,255,$transparency);
    @{$availableColors{ magenta }} = (255,0,255,$transparency);
    @{$availableColors{ cyan }} = (0,255,255,$transparency);
    #special 
    @{$availableColors{ orange }} = (255,165,0,$transparency);
    @{$availableColors{ blueGreen }} = (0,165,255,$transparency);
    #special 
    @{$availableColors{ yellowGreen }} = (0,255,165,$transparency);
    @{$availableColors{ redViolet }} = (255,0,165,$transparency);
    #special 
    @{$availableColors{ blueViolet }} = (165,0,255,$transparency);
    @{$availableColors{ yellowOrange }} = (165,255,0,$transparency);
    #special 
    @{$availableColors{ indigo }} = (75,0,130,$transparency);
    #special 
    @{$availableColors{ darkRed }} = (90,0,0,$transparency);
    @{$availableColors{ darkGreen }} = (0,90,0,$transparency);
    @{$availableColors{ darkBlue }} = (0,0,90,$transparency);
    #monoChromatic
    @{$availableColors{ monoRedStart }} = (255,203,203,$transparency);
    @{$availableColors{ monoRedEnd }} = (5,0,0,$transparency);
    @{$availableColors{ momoGreenStart }} = (233,255,203,$transparency);
    @{$availableColors{ monoGreenEnd }} = (0,5,0,$transparency);
    @{$availableColors{ monoBlueStart }} = (203,203,255,$transparency);
    @{$availableColors{ monoBlueEnd }} = (0,0,5,$transparency);
    
    #america 
    @{$availableColors{ lightRed }} = (255,165,165,$transparency);
    @{$availableColors{ lightBlue }} = (165,155,255,$transparency);
    
    # male/female
    @{$availableColors{ female }} = (255,0,255,$transparency);
    @{$availableColors{ male }} = (0,255,255,$transparency);
    
    # jon colors
    @{$availableColors{ jmb_orange }} = (255,105,0,$transparency);
    @{$availableColors{ jmb_darkOrange }} = (255,100,50,$transparency);
    @{$availableColors{ jmb_darkGreen }} = (25,100,50,$transparency);
    @{$availableColors{ jmb_gold }} = (255,215,0,$transparency);
    @{$availableColors{ jmb_lightRed }} = (255,128,128,$transparency);
    @{$availableColors{ jmb_darkRed }} = (128,0,0,$transparency);
    @{$availableColors{ jmb_lightBlue }} = (128,128,255,$transparency);
    @{$availableColors{ jmb_darkBlue }} = (0,0,128,$transparency);

    return(\%availableColors,$transparency);
}

=head2 initHeatmap

 Title     : initHeatmap
 Usage     : $img,$colorPalette,$nPosNegColorShades,$availableColors=initHeatmap(...)
 Function  : init a perl-GD image object
 Returns   : img object, colorPalette hash, nPosNegColorShades, availableColors hash
 Argument  : 

=cut

sub initHeatmap($$$;$$$) {
    # required
    my $imageHeight=shift;
    my $imageWidth=shift;
    my $colorString=shift;
    # optional
    my $missingColor="null";
    $missingColor=shift if @_;
    my $highlightColor="cyan";
    $highlightColor=shift if @_;
    my $verbose=0;
    $verbose=shift if @_;
    
    GD::Image->trueColor(1);
    my $img = new GD::Image($imageWidth,$imageHeight);
    
    #$img->saveAlpha(1);
    
    my $colorPalette={};
    my ($availableColors,$transparency)=getAvailableColors();
    my ($availableColorScales)=getAvailableColorScales();
    
    # first color becomes images background
    $colorPalette->{ BG } = $img->colorAllocate(200,200,200);
    
    $img->fill(0,0,$colorPalette->{ BG });
    
    # utility colors
    $colorPalette->{ W } = $img->colorAllocateAlpha(@{$availableColors->{ white }});
    $colorPalette->{ G } = $img->colorAllocateAlpha(@{$availableColors->{ grey }});
    $colorPalette->{ B } = $img->colorAllocateAlpha(@{$availableColors->{ black }});
    $colorPalette->{ N } = $img->colorAllocateAlpha(@{$availableColors->{ null }});
    
    # missing data color
    if($missingColor =~ /\./) {
        my @rgb = split(/\./, $missingColor);
        confess "color is not in the available color list nor is a valid rgb color ($missingColor)" if((@rgb != 3) and (@rgb != 4));
        foreach my $rgb_value (@rgb){
            confess "color is not in the available color list nor is a valid rgb color ($missingColor)" if(($rgb_value < 0) || ($rgb_value > 255));
        }
        my $tmpTransparency=$transparency;
        $tmpTransparency=$rgb[3] if(@rgb == 4);
        @{$availableColors->{ $missingColor }} = ($rgb[0],$rgb[1],$rgb[2],$tmpTransparency);
    } 
    confess "unknown color ($missingColor)" if(!(exists($availableColors->{ $missingColor })));
    $colorPalette->{ NA } = $img->colorAllocateAlpha(@{$availableColors->{ $missingColor }});
    
    # highlight data color
    if($highlightColor =~ /\./) {
        my @rgb = split(/\./, $highlightColor);
        confess "color is not in the available color list nor is a valid rgb color ($highlightColor)" if((@rgb != 3) and (@rgb != 4));
        foreach my $rgb_value (@rgb){
            confess "color is not in the available color list nor is a valid rgb color ($highlightColor)" if(($rgb_value < 0) || ($rgb_value > 255));
        }
        my $tmpTransparency=60;
        $tmpTransparency=$rgb[3] if(@rgb == 4);
        @{$availableColors->{ $highlightColor }} = ($rgb[0],$rgb[1],$rgb[2],$tmpTransparency);
    } else {
        confess "unknown color ($highlightColor)" if(!(exists($availableColors->{ $highlightColor })));
        my @rgb=@{$availableColors->{ $highlightColor }};
        $highlightColor=$highlightColor."_highlight";
        @{$availableColors->{ $highlightColor }} = ($rgb[0],$rgb[1],$rgb[2],60);
    }
    $colorPalette->{ HIGHLIGHT } = $img->colorAllocateAlpha(@{$availableColors->{ $highlightColor }});
    
    my $nColorShades = 255;
    $nColorShades -= keys(%{$colorPalette});
    
    my $nPositiveColorShades=floor($nColorShades/2);
    my $nNegativeColorShades=floor($nColorShades/2);
    
    my ($positiveColorString,$negativeColorString)=split(/___/,$colorString);
    
    # turn names color scale into names color list
    $positiveColorString=$availableColorScales->{ $positiveColorString } if(exists($availableColorScales->{ $positiveColorString }));
    $negativeColorString=$availableColorScales->{ $negativeColorString } if(exists($availableColorScales->{ $negativeColorString }));
    
    # positive color information
    $positiveColorString=$positiveColorString.",".$positiveColorString if($positiveColorString !~ /,/);
    my @positiveColorStringArray=split(/,/,$positiveColorString);
    my $nPositiveColors=@positiveColorStringArray;
    my $nPositiveShadesPerColor=floor($nPositiveColorShades/($nPositiveColors-1));    
    $nPositiveColorShades = (($nPositiveColors-1)*($nPositiveShadesPerColor-1));
    
    # negative color information
    $negativeColorString=$negativeColorString.",".$negativeColorString if($negativeColorString !~ /,/);
    my @negativeColorStringArray=split(/,/,$negativeColorString);
    my $nNegativeColors=@negativeColorStringArray;
    my $nNegativeShadesPerColor=floor($nNegativeColorShades/($nNegativeColors-1));    
    $nNegativeColorShades = (($nNegativeColors-1)*($nNegativeShadesPerColor-1));
    
    # adjust pos/neg color information to make symmetrical
    $nPositiveColorShades = $nNegativeColorShades if($nNegativeColorShades < $nPositiveColorShades);
    $nNegativeColorShades = $nPositiveColorShades if($nPositiveColorShades < $nNegativeColorShades);
    
    print STDERR "\tpositive color ($positiveColorString)\n" if($verbose);
    ($colorPalette,$img,$nPositiveColorShades)=calculateColorPalette($nPositiveColorShades,$positiveColorString,$availableColors,$colorPalette,"pc",$img,$transparency);

    print STDERR "\tnegative color ($negativeColorString)\n" if($verbose);
    ($colorPalette,$img,$nNegativeColorShades)=calculateColorPalette($nNegativeColorShades,$negativeColorString,$availableColors,$colorPalette,"nc",$img,$transparency);
        
    my $nPosNegColorShades=$nPositiveColorShades=$nNegativeColorShades;
    
    return($img,$colorPalette,$nPosNegColorShades,$availableColors);
}

=head2 calculateColorPalette

 Title     : calculateColorPalette
 Usage     : $colorPalette,$img,$nColorShades=calculateColorPalette(...)
 Function  : calculate a user-defined color palette [range]
 Returns   : nColorShades, colorString, availableColors hash, colorPalette hash, colorMode, img object, transparency value
 Argument  : 

=cut

sub calculateColorPalette($$$$$$$) {
    my $nColorShades=shift;
    my $colorString=shift;
    my $availableColors=shift;
    my $colorPalette=shift;
    my $colorMode=shift;
    my $img=shift;
    my $transparency=shift;
    
    my @colorStringArray=split(/,/,$colorString);
    my $nColors=@colorStringArray;
    confess "not enough colors! $nColors - $colorString" if($nColors < 2);

    my $shadesPerColor=floor($nColorShades/($nColors-1));    
    
    # JMB 08/14/14 change
    #validate that color either exists in palette or is a valid RGB color
    for(my $i=0;$i<$nColors;$i++) {
        my $currentColorName=$colorStringArray[$i];
        
        if($currentColorName =~ /\./) { # if color is RGB code
            my @rgb = split(/\./, $currentColorName);
            confess "color is not in the available color list nor is a valid rgb color ($currentColorName)" if((@rgb != 3) and (@rgb != 4));
            
            foreach my $rgb_value (@rgb){
                confess "color is not in the available color list nor is a valid rgb color ($currentColorName)" if(($rgb_value < 0) || ($rgb_value > 255));
            }
            
            # push user defined RGB code onto availableColors container
            my $tmpTransparency=$transparency;
            $tmpTransparency=$rgb[3] if(@rgb == 4);
            @{$availableColors->{ $currentColorName }} = ($rgb[0],$rgb[1],$rgb[2],$tmpTransparency);
            
        } else { # if color is named color
            confess "unknown color ($currentColorName)" if(!(exists($availableColors->{ $currentColorName })));
        }
    }
    

    for(my $i=0;$i<($nColors-1);$i++) {
        my $currentColorName=$colorStringArray[$i];
        my $currentColor=$availableColors->{ $currentColorName };
        my $currentColorRed=$currentColor->[0];
        my $currentColorGreen=$currentColor->[1];
        my $currentColorBlue=$currentColor->[2];
        
        my $nextColorName=$colorStringArray[$i+1];
        my $nextColor=$availableColors->{ $nextColorName };
        my $nextColorRed=$nextColor->[0];
        my $nextColorGreen=$nextColor->[1];
        my $nextColorBlue=$nextColor->[2];
        
        my $redDistance=($currentColorRed - $nextColorRed);
        my $greenDistance=($currentColorGreen - $nextColorGreen);
        my $blueDistance=($currentColorBlue - $nextColorBlue);        
        
        my $redColorStep=($redDistance/($shadesPerColor-1));
        my $greenColorStep=($greenDistance/($shadesPerColor-1));
        my $blueColorStep=($blueDistance/($shadesPerColor-1));
        
        for(my $c=0;$c<$shadesPerColor;$c++) {
            my $red=$currentColorRed;
            my $green=$currentColorGreen;
            my $blue=$currentColorBlue;
            
            $red -= round($c*$redColorStep);
            $green -= round($c*$greenColorStep);
            $blue -= round($c*$blueColorStep);
            
            $red=$nextColorRed if($c == ($shadesPerColor-1));
            $green=$nextColorGreen if($c == ($shadesPerColor-1));
            $blue=$nextColorBlue if($c == ($shadesPerColor-1));
            
            my $colorIndex=(($i*$shadesPerColor)+$c);
            $colorPalette->{ $colorMode }[$colorIndex] = $img->colorAllocateAlpha($red,$green,$blue,$transparency);
            
        }

    }
    
    $nColorShades = (($nColors-1)*$shadesPerColor);
    
    return($colorPalette,$img,$nColorShades);
}

=head2 getColorIndex

 Title     : getColorIndex
 Usage     : $colorIndex=getColorIndex(...)
 Function  : translate a signal to colorIndex
 Returns   : color index
 Argument  : score, colorScaleStart, colorScaleEnd, nColorShades, colorBucketSize

=cut

sub getColorIndex($$$$$) {
    my $score=shift;
    my $colorScaleStart=shift;
    my $colorScaleEnd=shift;
    my $nColorShades=shift;
    my $colorBucketSize=shift;
    
    $score=abs($score);
    
    my $colorIndex=0;
    if($score >= 0) { 
        if(($score >= $colorScaleStart) and ($score <= $colorScaleEnd)) {
            $score = $score - $colorScaleStart;
            $colorIndex = $score / $colorBucketSize if($colorBucketSize != 0);
            if($colorIndex > ($nColorShades-1)) { $colorIndex = ($nColorShades-1); } elsif($colorIndex < 0) { $colorIndex = 0; }
        } elsif($score < $colorScaleStart) {
            $colorIndex = 0;
        } elsif($score > $colorScaleEnd) {
            $colorIndex = ($nColorShades-1);
        } else {
            confess "error retreiving color index [$score] color: $colorScaleStart - $colorScaleEnd by $colorBucketSize ($nColorShades]";
        }
    }
    
    $colorIndex=floor($colorIndex);
    
    return($colorIndex);
}

=head2 chooseOptimalMissingValue

 Title     : chooseOptimalMissingValue
 Usage     : $num_zeros,$num_nans,$missingValue=chooseOptimalMissingValue(...)
 Function  : find optimal value for missing [sparse matrix format]
 Returns   : missing value
 Argument  : input matrix file path

=cut

sub chooseOptimalMissingValue($) {
    my $inputMatrix=shift;
    
    my $num_nans=0;
    my $num_zeros=0;
    
    my $lineNum=0;
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        if($lineNum > 0) { #header
            my @tmp=split(/\t/,$line);
            for(my $i=1;$i<@tmp;$i++) {
                if(($tmp[$i] eq "NA") or ($tmp[$i] eq "nan") or ($tmp[$i] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/))) {
                    $num_nans++;
                    next;
                } elsif($tmp[$i] == 0) {
                    $num_zeros++;
                    next;
                }
            }
        }
        $lineNum++;
   }
   close(IN);
   
    my $missingValue=0;
    if($num_nans > $num_zeros) {
        $missingValue = "NA";
    } else {
        $missingValue = 0;
    }
    
    return($num_zeros,$num_nans,$missingValue);

}

=head2 createTmpDir

 Title     : createTmpDir
 Usage     : $tmpDir=createTmpDir(...)
 Function  : create a tmp directory for tmp files used during processing
 Returns   : tmpDir file path
 Argument  : 

=cut

sub createTmpDir(;$) {
    #optional
    my $tmpDir="/tmp";
    $tmpDir=shift if @_;
    
    # remove trailing /
    $tmpDir =~ s/\/$//;
    
    my $uniq=getSmallUniqueString();
    $tmpDir = $tmpDir."/cWorld__".$uniq."/";
    
    system("mkdir -p '".$tmpDir."'");
    
    confess "could not create tmpDir (".$tmpDir.")!" if(!(-d($tmpDir)));
    
    return($tmpDir);
}

=head2 removeTmpDir

 Title     : removeTmpDir
 Usage     : removeTmpDir(...)
 Function  : garbage collect on a tmpDir
 Returns   : 
 Argument  : tmpDir file path

=cut

sub removeTmpDir($) {
    my $tmpDir=shift;
    
    confess "tmpDir does not exist! (".$tmpDir.")!" if(!(-d($tmpDir)));
    
    system("rm -rf '".$tmpDir."'") if($tmpDir =~ /\/cWorld__/);
}

=head2 headers2bed

 Title     : headers2bed
 Usage     : $headerBEDFile=headers2bed(...)
 Function  : generate a bed file of matrix headers
 Returns   : headerBED file path
 Argument  : matrixObject hash OR inc2header hash

=cut

sub headers2bed($) {
    my $input=shift;
    
    my ($matrixObject);
    if(-e $input) { # if this is a file
        $matrixObject=getMatrixObject($input);
    } elsif(exists($input->{ inc2header })) {
        $matrixObject=$input;
    } else {
        confess "invalid input, must be either matrixObject, or inputMatrix file";
    }
    
    my $inc2header=$matrixObject->{ inc2header };
    my $header2inc=$matrixObject->{ header2inc };
    my $numYHeaders=$matrixObject->{ numYHeaders };
    my $numXHeaders=$matrixObject->{ numXHeaders };
    my $numTotalHeaders=$matrixObject->{ numTotalHeaders };
    my $inputMatrixName=$matrixObject->{ inputMatrixName };
    my $NA_rowcols=$matrixObject->{ NArowcols };
    my $assembly=$matrixObject->{ assembly };
    
    my $headerBEDFile=$inputMatrixName."__".getSmallUniqueString().".headers.bed";
    open(BED,outputWrapper($headerBEDFile)) or confess "Could not open file [$headerBEDFile] - $!";
    
    print BED "track name='".$inputMatrixName."-headers' description='".$inputMatrixName."-headers' db=".$assembly." visibility=squish itemRgb=On useScore=0 color=0,0,0\n";
    
    my $enforceValidHeader=1;
    for(my $i=0;$i<$numTotalHeaders;$i++) {
        my $header=$inc2header->{ xy }->{$i};
        
        my $headerObject=getHeaderObject($header,$enforceValidHeader);
        
        my $headerChromosome="NA";
        $headerChromosome=$headerObject->{ chromosome } if(exists($headerObject->{ chromosome }));
        
        my $headerStart="NA";
        $headerStart=$headerObject->{ start } if(exists($headerObject->{ start }));
        my $headerEnd="NA";
        $headerEnd=$headerObject->{ end } if(exists($headerObject->{ end }));
        
        my $usableHeaderFlag=1;
        $usableHeaderFlag=0 if(exists($NA_rowcols->{$header}));
        
        # de group chromosome for UCSC use
        $headerChromosome=stripChromosomeGroup($headerChromosome);
        my $strand=".";
        
        my $color="0,0,0";
        $color="0,0,255" if($header =~ /_FOR_/);
        $color="255,0,0" if($header =~ /_REV_/);
        $color="0,255,255" if($header =~ /_LFOR_/);
        $color="255,0,255" if($header =~ /_LREV_/);
        $color="0,255,0" if($header =~ /_UNUSED_/);
        
        print BED "$headerChromosome\t$headerStart\t$headerEnd\t$header\t$usableHeaderFlag\t$strand\t$headerStart\t$headerEnd\t$color\n";
        
    }
    
    close(BED);
    
    confess "could not write BED file ($headerBEDFile)!" if(!(-e($headerBEDFile)));
    
    return($headerBEDFile);
    
}

=head2 combineBedFiles

 Title     : combineBedFiles
 Usage     : $bedFile=combineBedFiles(...)
 Function  : concatenate bed files, and collapse down to bed5
 Returns   : bed file path
 Argument  : array reference of bed file paths

=cut

sub combineBedFiles($;$) {
    # required
    my $bedFilesArrRef=shift;
    # optional
    my $bedName="";
    $bedName=shift if @_;
    
    $bedName = "__".$bedName if($bedName ne "");
    my $combinedBedFile=$bedName."__".getSmallUniqueString().".bed";
    
    open(OUT,outputWrapper($combinedBedFile)) or confess "Could not open file [$combinedBedFile] - $!";
    
    for(my $i=0;$i<@{$bedFilesArrRef};$i++) {
        my $bedFile = $bedFilesArrRef->[$i];
        
        open(IN,inputWrapper($bedFile)) or confess "Could not open file [$bedFile] - $!";
        while(my $line = <IN>) {
            chomp($line);
            next if(($line eq "") or ($line =~ m/^#/));
            
            # skip possible BED headers
            next if(($line eq "") or ($line =~ m/^#/));
            next if($line =~ /^#/);
            next if($line =~ /^track/);
            next if($line =~ /^chrom/);
            
            my @tmp=split(/\t/,$line);
            
            confess "ERROR-1: invalid BED format [$line]" if($line eq "");
            confess "ERROR-2: invalid BED format [$line]" if($tmp[0] !~ /^chr/);
            confess "ERROR-3: invalid BED format [$line]" if($tmp[1] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
            confess "ERROR-4: invalid BED format [$line]" if($tmp[2] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
            
            print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\n";
        }
        close(IN);
    }
    
    close(OUT);
    
    return($combinedBedFile);
}

=head2 midpointBedFile

 Title     : midpointBedFile
 Usage     : $bedFile=midpointBedFile(...)
 Function  : transform bed file into a bed file of the interval midpoints
 Returns   : bed file path
 Argument  : bed file path

=cut

sub midpointBedFile($$) {
    # required
    my $inputBedFile=shift;
    # optional
    my $bedName="";
    $bedName=shift if @_;
    
    $bedName = "__".$bedName if($bedName ne "");
    my $combinedBedFile=$bedName."__".getSmallUniqueString().".bed";
    
    my $numElements=0;
    
    open(OUT,outputWrapper($combinedBedFile)) or confess "Could not open file [$combinedBedFile] - $!";
    
    open(IN,inputWrapper($inputBedFile)) or confess "Could not open file [$inputBedFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        # skip possible BED headers
        next if(($line eq "") or ($line =~ m/^#/));
        next if($line =~ /^#/);
        next if($line =~ /^track/);
        next if($line =~ /^chrom/);
        
        my @tmp=split(/\t/,$line);
        
        confess "ERROR-1: invalid BED format [$line]" if($line eq "");
        confess "ERROR-2: invalid BED format [$line]" if($tmp[0] !~ /^chr/);
        confess "ERROR-3: invalid BED format [$line]" if($tmp[1] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        confess "ERROR-4: invalid BED format [$line]" if($tmp[2] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        
        my $midpoint=(($tmp[1]+$tmp[2])/2);
        
        $tmp[1]=floor($midpoint);
        $tmp[2]=ceil($midpoint);
        
        print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\n";
        $numElements++;
    }
    close(IN);
    close(OUT);
    
    return($combinedBedFile);
}

=head2 validateBED

 Title     : validateBED
 Usage     : validateBED(...)
 Function  : validate bed file is in correct UCSC format
 Returns   : live or die
 Argument  : bed file path

=cut

sub validateBED($) {
    my $bedFile=shift;
    
    # expect BED5+ format
    # chrom    chromStart    chromEnd    name    signalValue
    # track=
    
    # if pipe, do not process since reading pipe is destructive
    return if(-p $bedFile);
    
    open(IN,inputWrapper($bedFile)) or confess "Could not open file [$bedFile] - $!";
    while(my $line = <IN>) {
        chomp($line);        
        next if(($line eq "") or ($line =~ m/^#/));
        
        # skip possible BED headers
        next if(($line eq "") or ($line =~ m/^#/));
        next if($line =~ /^#/);
        next if($line =~ /^track/);
        next if($line =~ /^chrom/);
        
        my @tmp=split(/\t/,$line);
        
        confess "ERROR-1: invalid BED format [$line]" if($line eq "");
        confess "ERROR-2: invalid BED format [$line]" if($tmp[0] !~ /^chr/);
        confess "ERROR-3: invalid BED format [$line]" if($tmp[1] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        confess "ERROR-4: invalid BED format [$line]" if($tmp[2] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        
    }
    close(IN);
    
}

=head2 extendBED

 Title     : extendBED
 Usage     : $bedFile=extendBED(...)
 Function  : extend bed file intervals by X bp
 Returns   : bed file path
 Argument  : bed file path

=cut

sub extendBED($;$) {
    # required
    my $bedFile=shift;
    # optional
    my $elementExtension=0;
    $elementExtension=shift if @_;

    # expect BED5+ format
    # chrom    chromStart    chromEnd    name    signalValue
    # track=
    
    my $bedFileName=getFileName($bedFile);
    my $extendedBEDFile=$bedFileName.".exten".$elementExtension.".bed";
    
    open(OUT,outputWrapper($extendedBEDFile)) or confess "Could not open file [$extendedBEDFile] - $!";
    
    open(IN,inputWrapper($bedFile)) or confess "Could not open file [$bedFile] - $!";
    while(my $line = <IN>) {
        chomp($line);        
        next if(($line eq "") or ($line =~ m/^#/));
        
        # skip possible BED headers
        next if(($line eq "") or ($line =~ m/^#/));
        next if($line =~ /^#/);
        next if($line =~ /^track/);
        next if($line =~ /^chrom/);
        
        my @tmp=split(/\t/,$line);
        
        confess "ERROR-1: invalid BED format [$line]" if($line eq "");
        confess "ERROR-2: invalid BED format [$line]" if($tmp[0] !~ /^chr/);
        confess "ERROR-3: invalid BED format [$line]" if($tmp[1] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        confess "ERROR-4: invalid BED format [$line]" if($tmp[2] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        
        my $start=$tmp[1];
        my $end=$tmp[2];
        
        $start -= $elementExtension;
        $start=1 if($start < 1);
        $end += $elementExtension;
        
        $tmp[1]=$start;
        $tmp[2]=$end;
        
        print OUT join("\t", @tmp)."\n";
        
    }
    close(IN);
    
    close(OUT);
    
    return($extendedBEDFile);
    
}

=head2 intersectBED

 Title     : intersectBED
 Usage     : $bedFile=intersectBED(...)
 Function  : extend bed file intervals by X bp
 Returns   : bed file path
 Argument  : bed file path

=cut
    
sub intersectBED($$;$) {
    # required
    my $bedFile_1=shift; # header bed file
    my $bedFile_2=shift; # element bed file
    # optional
    my $elementExtension=0;
    $elementExtension=shift if @_;
    
    validateBED($bedFile_1);
    validateBED($bedFile_2);
    
    my $standardized_bedFile1=standardizeBED($bedFile_1);
    my $standardized_bedFile2=standardizeBED($bedFile_2);
    
    my $extended_bedFile2=extendBED($standardized_bedFile2,$elementExtension) if($elementExtension > 0);
    system("rm '".$standardized_bedFile2."'") if($elementExtension > 0);
    $standardized_bedFile2=$extended_bedFile2 if($elementExtension > 0);
    
    my $bedFileName_1=getFileName($standardized_bedFile1);
    my $bedFileName_2=getFileName($standardized_bedFile2);
    
    my $bedOverlapFile=$bedFileName_1."___".$bedFileName_2.".bed";
    
    system("bedtools intersect -a '".$standardized_bedFile1 ."' -b '".$standardized_bedFile2."' -wb > '".$bedOverlapFile."'");
    
    confess "could not write bed overlap file" if(!(-e($bedOverlapFile)));
    
    system("rm '".$standardized_bedFile1."'");
    system("rm '".$standardized_bedFile2."'");
    
    return($bedOverlapFile);
    
}

=head2 standardizeBED

 Title     : standardizeBED
 Usage     : $bedFile=standardizeBED(...)
 Function  : standardize bed file to bed 4 format
 Returns   : bed file path
 Argument  : bed file path

=cut

sub standardizeBED($) {
    # required
    my $bedFile=shift;
    
    # expect BED5+ format
    # chrom    chromStart    chromEnd    name    signalValue
    # track=
    
    my $bedFileName=getFileName($bedFile);
    $bedFileName=$bedFileName."__".getSmallUniqueString();
    my $standardizedBEDFile=$bedFileName.".bed4.bed";
    
    open(OUT,outputWrapper($standardizedBEDFile)) or confess "Could not open file [$standardizedBEDFile] - $!";
    
    my $lineNum=0;
    open(IN,inputWrapper($bedFile)) or confess "Could not open file [$bedFile] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        # skip possible BED headers
        next if($line =~ /^#/);
        next if($line =~ /^track/);
        next if($line =~ /^chrom/);
        
        my @tmp=split(/\t/,$line);
        
        confess "ERROR-1: invalid BED format [$line]" if($line eq "");
        confess "ERROR-2: invalid BED format [$line]" if($tmp[0] !~ /^chr/);
        confess "ERROR-3: invalid BED format [$line]" if($tmp[1] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        confess "ERROR-4: invalid BED format [$line]" if($tmp[2] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        
        $tmp[3] = $lineNum if(!defined($tmp[3]));
        print OUT "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\n";
        
        $lineNum++;
    }
    close(IN);
    
    close(OUT);
    
    return($standardizedBEDFile);
    
}

=head2 loadBED

 Title     : loadBED
 Usage     : $bedData=loadBED(...)
 Function  : standardize bed file to bed 4 format
 Returns   : bed file path
 Argument  : bed file path

=cut

sub loadBED($) {
    my $bedFile=shift;
    
    validateBED($bedFile);
    
    my @bedData=();
    my $bedInc=0;
    
    open(IN,inputWrapper($bedFile)) or confess "Could not open file [$bedFile] - $!";
    while(my $line = <IN>) {
        chomp($line);        
        next if(($line eq "") or ($line =~ m/^#/));
        
        # skip possible BED headers
        next if($line =~ /^#/);
        next if($line =~ /^track/);
        next if($line =~ /^chrom/);
        
        my @tmp=split(/\t/,$line);
        
        confess "ERROR-1: invalid BED format [$line]" if($line eq "");
        confess "ERROR-2: invalid BED format [$line]" if($tmp[0] !~ /^chr/);
        confess "ERROR-3: invalid BED format [$line]" if($tmp[1] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        confess "ERROR-4: invalid BED format [$line]" if($tmp[2] !~ (/^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/));
        
        my $coordinate=$tmp[0].":".$tmp[1]."-".$tmp[2];
        $bedData[$bedInc]{ chromosome }=$tmp[0];
        $bedData[$bedInc]{ start }=$tmp[1];
        $bedData[$bedInc]{ end }=$tmp[2];
        $bedData[$bedInc]{ name }=$tmp[3];
        
        my $name2="NA";
        $name2=$tmp[7] if(@tmp == 8);
        $bedData[$bedInc]{ name2 }=$name2;
        
        $bedData[$bedInc]{ coordinate }=$coordinate;
        $bedInc++;
        
    }
    close(IN);
    
    return(\@bedData);
}

=head2 transposeMatrix

 Title     : transposeMatrix
 Usage     : $matrixFile=transposeMatrix(...)
 Function  : transpose a matrix
 Returns   : input matrix file path
 Argument  : input matrix file path

=cut

sub transposeMatrix($;$) {
    # required
    my $inputMatrix=shift;
    # optional
    my $output=getSmallUniqueString().".transpose";
    $output = shift if @_;
    $output =~ s/\.matrix\.gz$//;
    $output =~ s/\.matrix$//;
    $output =~ s/\.gz$//;
    
    confess "inputMatrix [$inputMatrix] does not exist" if(!(-e $inputMatrix));

    my $paste_str="";
    my $file_str="";
    my $sub_paste_str="";
    my $sub_file_str="";

    my $tmpDir=createTmpDir();

    my $inputMatrixName=getFileName($inputMatrix);
    my $inputMatrixPath=getFilePath($inputMatrix);
    my $inputMatrixFile=$inputMatrixPath.$inputMatrixName;

    my $numLines=getNumberOfLines($inputMatrix);

    my $subNum=0;
    my $linenum=0;
    
    my $blockSize=50;
    
    open(IN,inputWrapper($inputMatrix)) or confess "Could not open file [$inputMatrix] - $!";
    while(my $line = <IN>) {
        chomp($line);
        next if(($line eq "") or ($line =~ m/^#/));
        
        #skip line is it is blank, starts with a # (means comment), starts with blank or starts with a character.
        next if(($line eq "") or ($line =~ m/^#/) or ($line =~ /^\s*$/));
        
        #turn all tabs into new lines
        $line =~ tr/\t/\n/;
        my $tmpfile=$tmpDir.$inputMatrixName.".".$linenum.".gz";
        
        #print the single row, which is now a column, into a temporary file.
        open(OUT,outputWrapper($tmpfile,0,1)) or confess "Could not open file [$tmpfile] - $!";
        print OUT $line;
        close(OUT);
        
        if($paste_str eq "") { $paste_str="<(gunzip -c '".$tmpfile."')"; } else { $paste_str=$paste_str." <(gunzip -c '".$tmpfile."')"; }
        if($file_str eq "") { $file_str="'".$tmpfile."'"; } else { $file_str=$file_str." '".$tmpfile."'"; }
        
        if(($linenum != 0) and (($linenum % $blockSize) == 0)) {
            my $subFile=$tmpDir.$inputMatrixName.".transpose.".$subNum.".gz";
            
            # run the paste
            my @args = ( "bash", "-c", "paste $paste_str | gzip > '$subFile'" );
            system(@args);
            # do the clean up
            system("rm $file_str");
            
            $subNum++;
            $paste_str="";
            $file_str="";
            
            if($sub_paste_str eq "") { $sub_paste_str="<(gunzip -c '".$subFile."')"; } else { $sub_paste_str=$sub_paste_str." <(gunzip -c '".$subFile."')"; }
            if($sub_file_str eq "") { $sub_file_str="'".$subFile."'"; } else { $sub_file_str=$sub_file_str." '".$subFile."'"; }
        }
        

        $linenum++;
    }
    close(IN);

    if($paste_str ne "") {
        my $subFile=$tmpDir.$inputMatrixName.".transpose.".$subNum.".gz";
        
        # run the paste
        my @args = ( "bash", "-c", "paste $paste_str | gzip > '$subFile'" );
        system(@args);
        # do the clean up
        system("rm $file_str");
        
        if($sub_paste_str eq "") { $sub_paste_str="<(gunzip -c '".$subFile."')"; } else { $sub_paste_str=$sub_paste_str." <(gunzip -c '".$subFile."')"; }
        if($sub_file_str eq "") { $sub_file_str="'".$subFile."'"; } else { $sub_file_str=$sub_file_str." '".$subFile."'"; }
    }

    my $transposeMatrix=$output.".matrix.gz";
    
    # run the paste
    my @args = ( "bash", "-c", "paste $sub_paste_str | gzip > '$transposeMatrix'" );
    system(@args);
    # do the clean up
    my @subFiles=split(/ /,$sub_file_str);
    for(my $f=0;$f<@subFiles;$f++) {
        my $sub_file=$subFiles[$f];
        system("rm '".$sub_file."'");
    }

    removeTmpDir($tmpDir);
    
    return($transposeMatrix);
}

=head2 flipBool

 Title     : flipBool
 Usage     : $bool=flipBool(...)
 Function  : flip a bool (0->1, 1->0)
 Returns   : bool
 Argument  : bool

=cut

sub flipBool($) {
    my $boolean=shift;
    
    confess "invalid bool value ($boolean)" if(($boolean != 0) and ($boolean != 1));
    
    return(1) if($boolean == 0);
    return(0) if($boolean);
}

=head2 outputWrapper

 Title     : outputWrapper
 Usage     : $outputFile=outputWrapper(...)
 Function  : wrap a outfile path to allow to inline gzip handing
 Returns   : outputfile string
 Argument  : outputfile string

=cut

sub outputWrapper($;$$) {
    # required
    my $outputFile=shift;
    # optional
    my $appendFlag=0;
    $appendFlag=shift if @_;
    my $suppressComments=0;
    $suppressComments=shift if @_;
    
    # disbale append flag if file does not yet exist
    $appendFlag = 0 if(!(-e $outputFile));
    
    my $outputCompressed=0;
    $outputCompressed=1 if($outputFile =~ /\.gz$/);

    my ($tmpOutputFile);
    
    if($appendFlag) {
        $tmpOutputFile = "| gzip -c >> '".$outputFile."'" if($outputCompressed);
        $tmpOutputFile = ">>".$outputFile if(!$outputCompressed);
    } else {
        $tmpOutputFile = "| gzip -c > '".$outputFile."'" if($outputCompressed);
        $tmpOutputFile = ">".$outputFile if(!$outputCompressed);
    }
    
    # disable comment(s)if (UCSC format file)
    $suppressComments = 1 if($outputFile =~ /\.bed$/);
    $suppressComments = 1 if($outputFile =~ /\.bedGraph$/);
    $suppressComments = 1 if($outputFile =~ /\.bed\.gz$/);
    $suppressComments = 1 if($outputFile =~ /\.bedGraph\.gz$/);
    $suppressComments = 1 if($outputFile =~ /\.wig$/);
    $suppressComments = 1 if($outputFile =~ /\.wig\.gz$/);
    
    if(!$suppressComments) {
        open(OUT,$tmpOutputFile) or confess "Could not open file [$tmpOutputFile] - $!";
        print OUT "## cworld::dekker\n";
        print OUT "## \n";
        print OUT "## Dekker Lab\n";
        print OUT "## Contact:\tBryan R. Lajoie\n";
        print OUT "## https://github.com/blajoie\n";
        print OUT "## \n";
        print OUT "## Version:\t".$cworld::dekker::VERSION."\n";
        print OUT "## Date:\t".getDate()."\n";
        print OUT "## Host:\t".getComputeResource()."\n";
        close(OUT);
    }
    
    $outputFile = "| gzip -c >> '".$outputFile."'" if($outputCompressed);
    $outputFile = ">>".$outputFile if(!$outputCompressed);
    
    return($outputFile);
}

=head2 inputWrapper

 Title     : inputWrapper
 Usage     : $inputfile=inputWrapper(...)
 Function  : wrap a infile path to allow to inline gzip handing
 Returns   : infile string
 Argument  : infile string

=cut

sub inputWrapper($) {
    my $inputFile=shift;
    
    $inputFile = "gunzip -c '".$inputFile."' | " if(($inputFile =~ /\.gz$/) and (!(-T($inputFile))));
    
    return($inputFile);
}

=head2 removeDiagonal

 Title     : removeDiagonal
 Usage     : $matrix=removeDiagonal(...)
 Function  : remove diagonal (x-y <= d), where d is num of diagonal to exclude.
 Returns   : matrix 2D hash
 Argument  : matrixObject hash, matrix 2D hash, excludeDiagonal value

=cut

sub removeDiagonal($$$) {
    my $matrixObject=shift;
    my $matrix=shift;
    my $excludeDiagonal=shift;
    
    my $inc2header=$matrixObject->{ inc2header };
    my $missingValue = $matrixObject->{ missingValue };
    
    my $numYHeaders=keys(%{$inc2header->{ y }});
    my $numXHeaders=keys(%{$inc2header->{ x }});

    for(my $y=0;$y<$numYHeaders;$y++) {
        for(my $x=0;$x<$numXHeaders;$x++) {
            
            if(abs($y-$x) <= $excludeDiagonal) {
                $matrix->{$y}->{$x}="NA" if($missingValue eq 0);
                delete($matrix->{$y}->{$x}) if($missingValue eq "NA");
            }
        }
    }
    
    return($matrix);
}

=head2 stitchMatrices

 Title     : stitchMatrices
 Usage     : $matrix=stitchMatrices(...)
 Function  : stitch together two symmetrical matrices, matrix 1 -> upper-diagonal, matrix 2 -> lower diagonal, diagonal -> NA
 Returns   : matrix 2D hash
 Argument  : matrixObject 1 hash, matrixObject 2 hash, matrix 1 2D hash, matrix 2 2D hash

=cut

sub stitchMatrices($$$$) {
    #required
    my $matrixObject_1=shift;
    my $matrixObject_2=shift;
    my $matrix_1=shift;
    my $matrix_2=shift;
    
    # assume matrix1 and matrix2 structure is same
    my $inc2header=$matrixObject_1->{ inc2header };
    
    my $numYHeaders=keys(%{$inc2header->{ y }});
    my $numXHeaders=keys(%{$inc2header->{ x }});

    my %stitchMatrix=();
    
    for(my $y=0;$y<$numYHeaders;$y++) {
        for(my $x=0;$x<$numXHeaders;$x++) {
            my $cScore="NA";
            if($y > $x) { #upper diagonal
                # use matrix_1 for upper diagonal
                $cScore=$matrixObject_1->{ missingValue };
                $cScore=$matrix_1->{$y}->{$x} if(defined($matrix_1->{$y}->{$x}));
            } elsif($x > $y) { # lower diagonal
                $cScore=$matrixObject_2->{ missingValue };
                $cScore=$matrix_2->{$y}->{$x} if(defined($matrix_2->{$y}->{$x}));
            } else { #exact diagonal
                $stitchMatrix{$y}{$x}="NA";
            }
            $stitchMatrix{$y}{$x}=$cScore if(($cScore ne "NA") and ($cScore ne "NA"));
        }
    }
    
    return(\%stitchMatrix);

}

=head2 intersectHeaders

 Title     : intersectHeaders
 Usage     : $elementHeaders=intersectHeaders(...)
 Function  : intersect two bed files
 Returns   : elementmatrix 2D hash
 Argument  : matrixObject 1 hash, matrixObject 2 hash, matrix 1 2D hash, matrix 2 2D hash

=cut

sub intersectHeaders($$;$) {
    # required
    my $matrixObject=shift;
    my $elementBedFile=shift;
    # optional
    my $elementExtension=0;
    $elementExtension=shift if @_;
    
    my $headerBEDFile=headers2bed($matrixObject);

    my $elementFileName=getFileName($elementBedFile);
    my $bedOverlapFile=intersectBED($headerBEDFile,$elementBedFile,$elementExtension);
    system("rm '".$headerBEDFile."'");

    
    validateBED($bedOverlapFile);
    
    my %elementHeaders=();      
    open(IN,inputWrapper($bedOverlapFile)) or confess "Could not open file [$bedOverlapFile] - $!";
    while(my $line = <IN>) {
        chomp($line);        
        next if(($line eq "") or ($line =~ m/^#/));
        
        # skip possible BED headers
        next if($line =~ /^#/);
        next if($line =~ /^track/);
        next if($line =~ /^chrom/);
        
        my @tmp=split(/\t/,$line);
        
        $elementHeaders{ $tmp[3] }=1;
    }
    close(IN);
    
    system("rm '".$bedOverlapFile."'");
    
    return(\%elementHeaders);
}

1; # End of cworld::dekker
