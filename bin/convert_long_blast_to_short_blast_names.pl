#!/usr/bin/perl -w

use strict;

while (<>) {
    next if /^#/; # mdb added 3/22/16 for LAST v731
    my @line = split(/\t/);
    my @item1 = split(/\|\|/, $line[0]);
    my @item2 = split(/\|\|/, $line[1]);
    $line[0] = $item1[6] if $item1[6];
    $line[1] = $item2[6] if $item2[6];
    print join ("\t", @line);
}
