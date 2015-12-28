# This file is a part of Julia. License is MIT: http://julialang.org/license

for (fmt, val) in (("%7.2f", "   1.23"),
                   ("%-7.2f", "1.23   "),
                   ("%07.2f", "0001.23"),
                   ("%.0f", "1"),
                   ("%#.0f", "1."),
                   ("%.4e", "1.2345e+00"))
    @test( @eval(@sprintf($fmt, big"1.2345") == $val))
end

# Inf / NaN handling
@test (@sprintf "%f" big"Inf") == "Inf"
@test (@sprintf "%f" big"NaN") == "NaN"

# scientific notation
@test (@sprintf "%.0e" big"3e142") == "3e+142"
@test (@sprintf "%#.0e" big"3e142") == "3.e+142"

for (val, res) in ((big"12345678.", "1.23457e+07"),
                   (big"1234567.8", "1.23457e+06"),
                   (big"123456.78", "123457"),
                   (big"12345.678", "12345.7"))
    @test( @sprintf("%.6g", val) == res)
end

for (fmt, val) in (("%10.5g", "     123.4"),
                   ("%+10.5g", "    +123.4"),
                   ("% 10.5g","     123.4"),
                   ("%#10.5g", "    123.40"),
                   ("%-10.5g", "123.4     "),
                   ("%-+10.5g", "+123.4    "),
                   ("%010.5g", "00000123.4"))
    @test( @eval(@sprintf($fmt, big"123.4") == $val))
end

@test( @sprintf( "%10.5g", big"-123.4" ) == "    -123.4")
@test( @sprintf( "%010.5g", big"-123.4" ) == "-0000123.4")
@test( @sprintf( "%.6g", big"12340000.0" ) == "1.234e+07")
@test( @sprintf( "%#.6g", big"12340000.0" ) == "1.23400e+07")
