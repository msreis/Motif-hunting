#!/usr/bin/perl -w

#------------------------------------------------------------------------------#
#                                                                              #
#  This program implements the IDENTIFY-SUBGRAPHS function described in        #
#  Figure 5 of the article:                                                    #
#                                                                              #
#  "A new approach for identification of cancer related pathways               #
#  using protein networks and genomic data"                                    #
#                                                                              #
#  André Fonseca, Marco Dimas Gubitoso, Marcelo S. Reis,                       #
#  Sandro José de Souza and Junior Barrera                                     #
#                                                                              #
#  Cancer Informatics 14(S5) 139–149 (2015).                                   #
#                                                                              #
#  doi: http://dx.doi.org/10.4137/CIN.S30800.                                  #
#                                                                              #
#  It also prints in the standard output the data depicted in Table 1.         #
#                                                                              #
#------------------------------------------------------------------------------#

#
#    Copyright (C) 2017 Marcelo S. Reis.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

use strict;
use Algorithm::Combinatorics qw(combinations);


# The Protein-Protein Interaction (PPI) graph depicted in Figure 4 of the paper. 
# When it is used as input for IDENTIFY-SUBGRAPHS, this function generates the 
# data showed in Table 1.
#
my %Pu = (
           7 => [6,5,3,2],
           6 => [],
           5 => [4],
           4 => [],
           3 => [],
           2 => [1],
           1 => [] 
         );

# Set X of vertices to be explored by the main algorithm (to reproduce the data
# showed in Table 1, it is necessary to constrain X to include vertex 7 only).
#
my %X = (7 => 1);

# Collection of sets of four vertices, in which the connected subgraphs of Pu 
# of size 4 will be stored.
#
my %G = ();

# Counter of the number of calls of the function.
#
my $call = 0;


# Implementation of IDENTIFY-SUBGRAPHS. 
#
# It receives: - a PPI graph Pu;
#              - a subset X of the vertices of Pu; 
#              - a vertex v in X;
#              - an integer k, 1 <= k < 4.
#
# It returns:  - all connected subgraphs of Pu that are induced by the vertices 
#                in a union of X and a set S of k vertices such that: 
#                  * for each vertex s in S, S < v;
#                  * for each vertex s in S, there is a path between s and v;
#                  * |X \cup S| = 4.
#
sub IDENTIFY_SUBGRAPHS
{
  my ($Pu, $G, $X, $v, $k) = @_;

  $call++;
  printf "%7s  ,  %55s  ,  %7s  ,  %d  ,  %d\n", 
    $call, print_collection ($G), 
    "{" . join(",", reverse sort keys %{$X}) . "}",
    $v, $k;
    
  if (scalar @{$Pu->{$v}} >= $k) # if the number of neighbors of $v is smaller
                                 # than $k, then the number of sets $S taken
                                 # $k at a time is zero.
  {
    my $iterator = combinations($Pu->{$v}, $k);
    while (my $S = $iterator->next)
    {
      my %S_cup_X = map { $_ => 1 } @{$S};      
      $S_cup_X{$_} = 1 foreach keys %{$X};      
      
      if (scalar keys %S_cup_X == 4)
      {
        $G{"{" . join(",", reverse sort keys %S_cup_X) . "}"} = 1;
      }
      else
      {      
        foreach my $s (@{$S})
        {
          if (! defined $X->{$s})
          {
            IDENTIFY_SUBGRAPHS 
              ($Pu, $G, \%S_cup_X, $s, 4 - scalar keys %S_cup_X);
          }
        }
      }
    } # end while
  } # end if
  
  if ($k - 1 > 0)
  {
    IDENTIFY_SUBGRAPHS ($Pu, $G, $X, $v, $k - 1); 
  }

} # end sub


# Subroutine that prints a collection G
#
sub print_collection
{
  my $G  = $_[0];
  my $string = "{";
  if (scalar keys %{$G} > 0)
  {
    foreach my $set (reverse sort keys %{$G})
    {
      $string .=  $set . ", ";
    }
    $string = substr ($string, 0, -2);
  }
  return $string . "}";
}


# Main program (it makes the initial call of the function).
#
printf " # call  ,  %55s  , " .
      " %7s  ,  v  ,  k\n", 
      "G (only the vertex sets of the graphs are shown)", "X";

IDENTIFY_SUBGRAPHS (\%Pu, \%G, \%X, 7, 3);

printf "%7s  ,  %55s  ,  %7s  ,  %s  ,  %s\n", 
  "**", print_collection (\%G), "-", "-", "-";

# End of program.
#
exit 0;

