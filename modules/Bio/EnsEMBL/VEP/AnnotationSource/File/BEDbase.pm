=head1 LICENSE

Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

# EnsEMBL module for Bio::EnsEMBL::VEP::AnnotationSource::File::BED
#
#

=head1 NAME

Bio::EnsEMBL::VEP::AnnotationSource::File::BEDbase - annotation source for BEDs with delimited per-base scores

=head1 SYNOPSIS

my $as = Bio::EnsEMBL::VEP::AnnotationSource::File::BEDbase->new({
  config => $config,
  file   => "my_features.bed.gz",
  type   => "exact"
});

$as->annotate_InputBuffer($ib);

=head1 DESCRIPTION

Modified BED format custom annotation source. BED files must be chromosome/pos
sorted, compressed with bgzip and indexed with tabix.

Each line of BED represents a chromosomal range, with the 4th column containing
scores or values separated by commas, one value per base in the range start-end.

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::VEP::AnnotationSource::File::BEDbase;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Variation::Utils::VariationEffect qw(overlap);
use Bio::EnsEMBL::IO::Parser::BedTabix;

use base qw(Bio::EnsEMBL::VEP::AnnotationSource::File::BED);


=head2 _create_records
 
  Example    : $records = $as->_create_records();
  Description: Create a custom annotation record from the current
               record as read from the annotation source.
  Returntype : arrayref of hashrefs
  Exceptions : none
  Caller     : annotate_VariationFeature()
  Status     : Stable

=cut

sub _create_records {
  return [{name => join(",", @{$_[1]})}];
}


=head2 _record_overlaps_VF
 
  Arg 1      : Bio::EnsEMBL::Variation::VariationFeature
  Example    : $overlap_ok = $as->_record_overlaps_VF($vf);
  Description: Determine whether the given VariationFeature overlaps
               the current record, depending on the set type()
  Returntype : bool
  Exceptions : none
  Caller     : annotate_VariationFeature()
  Status     : Stable

=cut

sub _record_overlaps_VF {
  my $self = shift;
  my $vf = shift;

  my $parser = $self->parser();
  my ($vs, $ve) = ($vf->{start}, $vf->{end});

  $DB::single = 1;

  # can't do insertions
  return 0 unless $vs <= $ve;

  my $ps = $parser->get_start;

  if(overlap($ps, $parser->get_end, $vs, $ve)) {
    my @values = split(",", $parser->get_name);
    return [map {$values[$_+1]} grep {$_ >= 0 && $_ <= $#values} (($vs - $ps)..($ve - $ps))]; 
  }
  else {
    return 0;
  }
}

1;