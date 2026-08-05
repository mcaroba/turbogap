#!/usr/bin/env bash

# Exit on error
set -e

usage() {
  echo "Usage: $0 [-i] <file_or_directory>"
  echo "  -i  In-place editing (overwrites original file)"
  echo "  -h  Show this help message"
  exit 1
}

IN_PLACE=0

while getopts "ih" opt; do
  case ${opt} in
  i) IN_PLACE=1 ;;
  h) usage ;;
  *) usage ;;
  esac
done
shift $((OPTIND - 1))

TARGET="$1"

if [ -z "$TARGET" ]; then
  echo "Error: Missing target file or directory."
  usage
fi

process_file() {
  local file="$1"

  local perl_script='
    use strict;
    use warnings;

    local $/ = undef;
    my $content = <STDIN>;

    # Split file into line arrays to process sequential blocks without overlapping multiline regex
    my @lines = split /\r?\n/, $content;
    my @output = ();
    my $i = 0;

    while ($i < scalar @lines) {
        my $line = $lines[$i];

        # Detect Fortran declaration containing "::"
        if ($line =~ /^([ \t]*[a-zA-Z0-9_(),\s*&:]+::[ \t]*)(.*)/) {
            my $hdr = $1;
            my $body = $2;

            # 1. Gather all continuation lines using &
            while ($body =~ /&[ \t]*$/ || ($i + 1 < scalar @lines && $lines[$i+1] =~ /^[ \t]*&/)) {
                $body =~ s/&[ \t]*$//;
                $i++;
                last if $i >= scalar @lines;
                my $next = $lines[$i];
                $next =~ s/^[ \t]*&[ \t]*//;
                $body .= " " . $next;
            }

            # 2. Clean header: remove stray & and fix duplicated keywords
            $hdr =~ s/&//g;
            $hdr =~ s/\s+/ /g;
            while ($hdr =~ s/\b(\w+)\s*::\s*\1\b/$1/gi) {}
            $hdr =~ s/(::\s*)+/:: /g;

            # 3. Clean body whitespace
            $body =~ s/^\s+|\s+$//g;

            # 4. Safely split variable list by comma (respecting nested parentheses)
            my @vars = ();
            my $depth = 0;
            my $cur = "";
            for my $char (split //, $body) {
                if ($char eq "(") { $depth++; }
                elsif ($char eq ")") { $depth--; }

                if ($char eq "," && $depth == 0) {
                    $cur =~ s/^\s+|\s+$//g;
                    push @vars, $cur if $cur ne "";
                    $cur = "";
                } else {
                    $cur .= $char;
                }
            }
            $cur =~ s/^\s+|\s+$//g;
            push @vars, $cur if $cur ne "";

            # 5. Push formatted individual declarations
            for my $var (@vars) {
                push @output, $hdr . $var;
            }
        } else {
            push @output, $line;
        }
        $i++;
    }

    print join("\n", @output) . "\n";
    '

  if [ "$IN_PLACE" -eq 1 ]; then
    local temp_file="${file}.tmp"
    perl -e "$perl_script" <"$file" >"$temp_file"
    mv "$temp_file" "$file"
    echo "Processed (in-place): $file"
  else
    perl -e "$perl_script" <"$file"
  fi
}

if [ -f "$TARGET" ]; then
  process_file "$TARGET"
elif [ -d "$TARGET" ]; then
  find "$TARGET" -type f \( -name "*.f90" -o -name "*.F90" -o -name "*.f" -o -name "*.f03" \) | while read -r f90_file; do
    process_file "$f90_file"
  done
else
  echo "Error: '$TARGET' is not a valid file or directory."
  exit 1
fi
