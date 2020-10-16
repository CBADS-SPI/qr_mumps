my $intmpl = 0;
my $tmpl   = "";
my $key    = "";
my %templates = ();

while (<>) {
    my $line = $_; # or simply "print;"

    if ($intmpl eq 0) {
        if (index($line, "!beg_template:") != -1) {
            my @words = split / /, $line;
            $key = @words[1];
            $key =~ s/\R//g;
            # print "Found template for ",@words[1];
            $intmpl = 1;
            $tmpl = "";
        } 
    } else {
        if (index($line, "!end_template") != -1) {
            $intmpl = 0;
            # print "\n\nHere is a template for ",$key,"\n";
            # print $tmpl;
            $templates{$key} .= $tmpl;
        } else {
            $tmpl = $tmpl.$line;
        }
    } 
}



@arrays = (
    { pref => "p", check => "associated", array => "pointer"    , move => "dst => src; nullify(src)"  },
    { pref => "a", check => "allocated",  array => "allocatable", move => "call move_alloc(src, dst)" }
    );

@ariths = (
    { arith => "s" , type => "real",    prec => "(kind(1.e0))" },
    { arith => "d" , type => "real",    prec => "(kind(1.d0))" },
    { arith => "c" , type => "complex", prec => "(kind(1.e0))" },
    { arith => "z" , type => "complex", prec => "(kind(1.d0))" },
    { arith => "i" , type => "integer", prec => "(kind=4)"     },
    { arith => "i8", type => "integer", prec => "(kind=8)"     }
    );

@sizes = (
    { rank => "1", dim => ":"    , sizes => "m"   ,    mult => "m"    , checksize => "m .le. 0" },
    { rank => "2", dim => ":,:"  , sizes => "m, n",    mult => "m*n"  , checksize => "(m .le. 0) .or. (n .le. 0)" },
    { rank => "3", dim => ":,:,:", sizes => "m, n, k", mult => "m*n*k", checksize => "(m .le. 0) .or. (n .le. 0) .or. (k .le. 0)" }
    );

# my @infaces = ("alloc", "dealloc", "move_alloc", "size", "realloc" );
# my %templates = (alloc => $talloc, dealloc => $tdealloc, move_alloc => $tmovealloc, size => $tsize, realloc => $trealloc);

print "module qrm_mem_mod\n"  ;
print "  use iso_c_binding\n"  ;
print "  use qrm_error_mod\n"  ;
print "  use qrm_const_mod\n"  ;
print "  use qrm_memhandling_mod\n";
print "  use qrm_pthread_mod\n\n";


my $ns = scalar(@sizes);
    
foreach my $key ( keys %templates )
{
    print "  interface qrm_",$key,"\n";
    foreach my $arr (@arrays)
    {
        foreach my $ari (@ariths)
        {
            print "    module procedure ";
            my $cnt = 0;
            foreach my $size (@sizes)
            {
                my %hash = (%{$arr}, %{$ari}, %{$size});
                print " qrm_",$hash{"pref"},$key,"_",$hash{"rank"},$hash{"arith"};
                if($key eq "realloc")
                {
                    last;
                }
                if($cnt < $ns-1)
                {
                    print ",";
                }
                $cnt = $cnt+1;
            }
            print "\n";
        }
    }
    print "  end interface\n\n";
}




print "\n\ncontains\n\n";


foreach my $key ( keys %templates )
{
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"; 
    print "!                     ",$key," routines                           !\n"; 
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"; 
    foreach my $arr (@arrays)
    {
        foreach my $ari (@ariths)
        {
            foreach my $size (@sizes)
            {
                my %hash = (%{$arr}, %{$ari}, %{$size});
                # print Dumper(\%hash), "\n";
                $template = $templates{$key};
                $template =~ s/\#(.*?)\#/$hash{"\L$1"}/g;
                print $template, "\n";
                if($key eq "realloc")
                {
                    last;
                }
            }
        }
    }
}

print "\nend module qrm_mem_mod"
