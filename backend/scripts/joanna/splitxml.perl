#!/usr/bin/env perl

cpam XML::Twig
download file from 

# - - - - - - - - -


use strict;
use warnings;

use XML::Twig;

#new document. Manually set xmlns - could copy this from 'original'
#instead though.
my $new_doc = XML::Twig->new;
$new_doc->set_root(
   XML::Twig::Elt->new(
      'urlset', { xmlns => "http://www.hmdb.ca" }
   )
);
$new_doc->set_pretty_print('indented_a');

my $elt_count    = 0;
my $elts_per_doc = 2;
my $count_of_xml = 0;

#handle each 'url' element.
sub handle_url {
   my ( $twig, $elt ) = @_;
   #more than the count, we output this doc, close it,
   #then create a new one.
   if ( $elt_count >= $elts_per_doc ) {
      $elt_count = 0;
      open( my $output, '>', "new_xml_" . $count_of_xml++ . ".xml" )
        or warn $!;
      print {$output} $new_doc->sprint;
      close($output);
      $new_doc = XML::Twig->new();
      $new_doc->set_root(
         XML::Twig::Elt->new(
            'urlset',
            { xmlns => "http://www.hmdb.ca" }
         )
      );
      $new_doc->set_pretty_print('indented_a');
   }
   #cut this element, paste it into new doc.
   #note - this doesn't alter the original on disk - only the 'in memory'
   #copy.
   $elt->cut;
   $elt->paste( $new_doc->root );
   $elt_count++;
   #purge clears any _closed_ tags from memory, so it preserves
   #structure.
   $twig->purge;
}

#set a handler, start the parse.

my $twig = XML::Twig->new( twig_handlers => { 'url' => \&handle_url } ) ->parsefile ( 'your_file.xml' );
